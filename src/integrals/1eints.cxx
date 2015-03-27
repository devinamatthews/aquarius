#include "1eints.hpp"

#include "kei.hpp"
#include "ovi.hpp"
#include "nai.hpp"

/**
 * Compute the index of a function in cartesian angular momentum.
 */
#define FUNC_CART(x,y,z) ((((x)*(3+(x)+2*((y)+(z))))/2) + (y))
/**
 * Compute the index of a function in spherical harmonic angular momentum.
 *
 * Regular spherical harmonics are referenced by n=l, l>=m>=-l. Contaminants may also be referenced by
 * n>l>=0, n-l even.
 */
#define FUNC_SPHER(n,l,m) ((((n)-(l))*((n)+(l)-1))/2 + 2*(n) + ((m) > 0 ? -2*(m) : 2*(m)+1))

#define IDX_EQ(i,r,e,j,s,f) ((i) == (j) && (r) == (s) && (e) == (f))
#define IDX_GE(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >= (f)))))
#define IDX_GT(i,r,e,j,s,f) ((i) > (j) || ((i) == (j) && ((r) > (s) || ((r) == (s) && (e) >  (f)))))

using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;
using namespace aquarius::tensor;

namespace aquarius
{
namespace integrals
{

OneElectronIntegrals::OneElectronIntegrals(const Shell& a, const Shell& b)
: sa(a), sb(b), group(a.getCenter().getPointGroup()),
  ca(a.getCenter()), cb(b.getCenter()), la(a.getL()), lb(b.getL()),
  na(a.getNPrim()), nb(b.getNPrim()), ma(a.getNContr()), mb(b.getNContr()),
  da(a.getDegeneracy()), db(b.getDegeneracy()), fsa(a.getNFunc()), fsb(b.getNFunc()),
  za(a.getExponents()), zb(b.getExponents()), num_processed(0)
{
    fca = (la+1)*(la+2)/2;
    fcb = (lb+1)*(lb+2)/2;

    int nints = 0;
    for (int j = 0;j < fsb;j++)
    {
        for (int i = 0;i < fsa;i++)
        {
            for (int s = 0;s < db;s++)
            {
                for (int r = 0;r < da;r++)
                {
                    const Representation& w = group.getIrrep(a.getIrrepOfFunc(i, r));
                    const Representation& x = group.getIrrep(b.getIrrepOfFunc(j, s));

                    if ((w*x).isTotallySymmetric()) nints += ma*mb;
                }
            }
        }
    }

    ints.resize(nints);
}

void OneElectronIntegrals::run()
{
    so(ints.data());
}

size_t OneElectronIntegrals::process(const Context& ctx, const vector<int>& idxa, const vector<int>& idxb,
                                     size_t nprocess, double* integrals, idx2_t* indices, double cutoff)
{
    const PointGroup& group = ca.getPointGroup();

    size_t m = 0;
    size_t n = 0;
    for (int j = 0;j < fsb;j++)
    {
        for (int i = 0;i < fsa;i++)
        {
            for (int s = 0;s < db;s++)
            {
                for (int r = 0;r < da;r++)
                {
                    const Representation w = group.getIrrep(sa.getIrrepOfFunc(i,r));
                    const Representation x = group.getIrrep(sb.getIrrepOfFunc(j,s));

                    if (!(w*x).isTotallySymmetric()) continue;

                    for (int f = 0;f < mb;f++)
                    {
                        for (int e = 0;e < ma;e++)
                        {
                            if (num_processed > m)
                            {
                                m++;
                                continue;
                            }

                            if (aquarius::abs(ints[m]) > cutoff && (&sa != &sb || IDX_GE(i,r,e,j,s,f)))
                            {
                                indices[n].i = sa.getIndex(ctx, idxa, i, e, r);
                                indices[n].j = sb.getIndex(ctx, idxb, j, f, s);
                                integrals[n++] = ints[m];
                            }

                            num_processed++;
                            m++;

                            if (n >= nprocess) return n;
                        }
                    }
                }
            }
        }
    }

    return n;
}

void OneElectronIntegrals::prim(const vec3& posa, int e,
                                const vec3& posb, int f, double* integrals)
{
    assert(0);
}

void OneElectronIntegrals::prims(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    #pragma omp parallel
    {
        #pragma omp for
        for (int m = 0;m < na*nb;m++)
        {
            int f = m/na;
            int e = m%na;
            prim(posa, e, posb, f, integrals+fca*fcb*m);
        }
    }
}

void OneElectronIntegrals::contr(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    vector<double> pintegrals(fca*fcb*na*nb);
    prims(posa, posb, pintegrals.data());
    prim2contr2r(fca*fcb, pintegrals.data(), integrals);
}

void OneElectronIntegrals::spher(const vec3& posa, const vec3& posb,
                                 double* integrals)
{
    vector<double> cintegrals(fca*fcb*na*nb);
    contr(posa, posb, integrals);
    cart2spher2r(ma*mb, integrals, cintegrals.data());
    transpose(fsa*fsb, ma*mb, 1.0, cintegrals.data(), fsa*fsb,
                              0.0,  integrals       ,   ma*mb);
}

void OneElectronIntegrals::so(double* integrals)
{
    const PointGroup& group = ca.getPointGroup();

    vector<double> aointegrals(fca*fcb*na*nb);

    int lambdar;
    vector<int> dcrr = group.DCR(ca.getStabilizer(), cb.getStabilizer(), lambdar);
    double coef = (double)group.getOrder()/(double)lambdar;

    for (int r : dcrr)
    {
        spher(ca.getCenter(0), cb.getCenter(cb.getCenterAfterOp(r)), aointegrals.data());
        scal(aointegrals.size(), coef, aointegrals.data(), 1);
        ao2so2(ma*mb, r, aointegrals.data(), integrals);
    }
}

void OneElectronIntegrals::ao2so2(size_t nother, int r, double* aointegrals, double* sointegrals)
{
    for (int j = 0;j < fsb;j++)
    {
        for (int i = 0;i < fsa;i++)
        {
            for (int f = 0;f < db;f++)
            {
                for (int e = 0;e < da;e++)
                {
                    int w = sa.getIrrepOfFunc(i,e);
                    int x = sb.getIrrepOfFunc(j,f);

                    if (!(group.getIrrep(w)*group.getIrrep(x)).isTotallySymmetric()) continue;

                    double fac = sb.getParity(j,r)*group.character(x,r);
                    axpy(nother, fac, aointegrals, 1, sointegrals, 1);

                    sointegrals += nother;
                }
            }
            aointegrals += nother;
        }
    }
}

void OneElectronIntegrals::cart2spher2r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    m = fsb;
    n = fca*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,j]' x [xa,b]' = [j,xa]
        gemm('T', 'T', m, n, k, 1.0, sb.getCart2Spher().data(), k, buf1, n, 0.0, buf2, m);
    }
    else
    {
        // [xa,j]' = [j,xa]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    m = fsa;
    n = fsb*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,i]' x [jx,a]' = [i,jx]
        gemm('T', 'T', m, n, k, 1.0, sa.getCart2Spher().data(), k, buf2, n, 0.0, buf1, m);
    }
    else
    {
        // [jx,i]' = [i,jx]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::cart2spher2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    n = fsa;
    m = fcb*nother;
    k = fca;

    if (sa.isSpherical())
    {
        // [a,bx]' x [a,i]  = [bx,i]
        gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCart2Spher().data(), k, 0.0, buf2, m);
    }
    else
    {
        // [i,bx]' = [bx,i]
        transpose(n, m, 1.0, buf1, n, 0.0, buf2, m);
    }

    n = fsb;
    m = fsa*nother;
    k = fcb;

    if (sb.isSpherical())
    {
        // [b,xi]' x [b,j] = [xi,j]
        gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCart2Spher().data(), k, 0.0, buf1, m);
    }
    else
    {
        // [j,xi]' = [xi,j]
        transpose(n, m, 1.0, buf2, n, 0.0, buf1, m);
    }

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2r(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [b,j]' x [xa,b]' = [j,xa]
    m = mb;
    n = na*nother;
    k = nb;
    gemm('T', 'T', m, n, k, 1.0, sb.getCoefficients().data(), k, buf1, n, 0.0, buf2, m);

    // [a,i]' x [jx,a]' = [i,jx]
    m = ma;
    n = mb*nother;
    k = na;
    gemm('T', 'T', m, n, k, 1.0, sa.getCoefficients().data(), k, buf2, n, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

void OneElectronIntegrals::prim2contr2l(size_t nother, double* buf1, double* buf2)
{
    size_t m, n, k;

    // [a,bx]' x [a,i] = [bx,i]
    m = ma;
    n = nb*nother;
    k = na;
    gemm('T', 'N', m, n, k, 1.0, buf1, k, sb.getCoefficients().data(), k, 0.0, buf2, m);

    // [b,xi]' x [b,j] = [xi,j]
    m = mb;
    n = ma*nother;
    k = nb;
    gemm('T', 'N', m, n, k, 1.0, buf2, k, sa.getCoefficients().data(), k, 0.0, buf1, m);

    copy(m*n, buf1, 1, buf2, 1);
}

OneElectronIntegralsTask::OneElectronIntegralsTask(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    addProduct(Product("ovi", "S", reqs));
    addProduct(Product("kei", "T", reqs));
    addProduct(Product("nai", "G", reqs));
    addProduct(Product("1ehamiltonian", "H", reqs));
}

bool OneElectronIntegralsTask::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");

    Context ctx(Context::ISCF);

    const vector<int>& N = molecule.getNumOrbitals();
    int n = molecule.getGroup().getNumIrreps();

    vector<int> irrep;
    for (int i = 0;i < n;i++) irrep += vector<int>(N[i],i);

    vector<uint16_t> start(n,0);
    for (int i = 1;i < n;i++) start[i] = start[i-1]+N[i-1];

    vector<vector<int> > idx = Shell::setupIndices(ctx, molecule);
    vector<Shell> shells(molecule.getShellsBegin(), molecule.getShellsEnd());
    vector<vector<tkv_pair<double> > > ovi_pairs(n), nai_pairs(n), kei_pairs(n);
    vector<Center> centers;

    for (auto& atom : molecule.getAtoms())
    {
        centers.push_back(atom.getCenter());
    }

    int block = 0;
    for (int a = 0;a < shells.size();++a)
    {
        for (int b = 0;b <= a;++b)
        {
            if (block%arena.size == arena.rank)
            {
                IshidaOVI s(shells[a], shells[b]);
                IshidaKEI t(shells[a], shells[b]);
                IshidaNAI g(shells[a], shells[b], centers);

                s.run();
                t.run();
                g.run();

                size_t nint = s.getIntegrals().size();
                vector<double> ints(nint);
                vector<idx2_t> idxs(nint);
                size_t nproc;

                nproc = s.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];
    
                                ovi_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) ovi_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }

                nproc = t.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];

                                kei_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) kei_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }

                nproc = g.process(ctx, idx[a], idx[b], nint, ints.data(), idxs.data());
                for (int k = 0;k < nproc;k++)
                {
                    int irr = irrep[idxs[k].i];
                    assert(irr == irrep[idxs[k].j]);

                    uint16_t i = idxs[k].i-start[irr];
                    uint16_t j = idxs[k].j-start[irr];

                                nai_pairs[irr].push_back(tkv_pair<double>(i*N[irr]+j, ints[k]));
                    if (i != j) nai_pairs[irr].push_back(tkv_pair<double>(j*N[irr]+i, ints[k]));
                }
            }

            block++;
        }
    }

    OVI *ovi = new OVI(arena, molecule.getGroup(), N);
    KEI *kei = new KEI(arena, molecule.getGroup(), N);
    NAI *nai = new NAI(arena, molecule.getGroup(), N);
    OneElectronHamiltonian *oeh = new OneElectronHamiltonian(arena, molecule.getGroup(), N);

    for (int i = 0;i < n;i++)
    {
        vector<int> irreps(2,i);
        (*ovi).writeRemoteData(irreps, ovi_pairs[i]);
        (*kei).writeRemoteData(irreps, kei_pairs[i]);
        (*nai).writeRemoteData(irreps, nai_pairs[i]);
        (*oeh).writeRemoteData(irreps, kei_pairs[i]);
        (*oeh).writeRemoteData(irreps, 1.0, 1.0, nai_pairs[i]);
    }

    put("S", ovi);
    put("T", kei);
    put("G", nai);
    put("H", oeh);

    return true;
}

}
}

REGISTER_TASK(aquarius::integrals::OneElectronIntegralsTask,"1eints");
