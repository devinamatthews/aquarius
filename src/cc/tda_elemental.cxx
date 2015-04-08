#include "tda_elemental.hpp"

using namespace El;

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace cc
{

template <typename U>
ElementalTDA<U>::ElementalTDA(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("tda.TDAevals", "TDAevals", reqs));
    addProduct(Product("tda.TDAevecs", "TDAevecs", reqs));
}

template <typename U>
bool ElementalTDA<U>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    int nirrep = group.getNumIrreps();

    const auto& W = get<TwoElectronOperator<U>>("H");
    const Space& occ = W.occ;
    const Space& vrt = W.vrt;

    SpinorbitalTensor<U> Hguess("Hguess", arena, group, {vrt,occ}, {1,1}, {1,1});
    Hguess = 0;

    const SpinorbitalTensor<U>& FAB = W.getAB();
    const SpinorbitalTensor<U>& FIJ = W.getIJ();
    const SpinorbitalTensor<U>& WAIBJ = W.getAIBJ();

    Hguess["aibi"]  = FAB["ab"];
    Hguess["aibj"] -= WAIBJ["aibj"];
    Hguess["aiaj"] -= FIJ["ij"];

    auto& TDAevecs = put("TDAevecs", new vector<unique_vector<SpinorbitalTensor<U>>>(nirrep));
    auto& TDAevals = put("TDAevals", new vector<vector<U>>(nirrep));

    for (int R = 0;R < nirrep;R++)
    {
        const Representation& irr_R = group.getIrrep(R);

        int ntot = 0;
        for (int i = 0, count = 0;i < nirrep;i++)
        {
            const Representation& irr_i = group.getIrrep(i);
            for (int a = 0;a < nirrep;a++)
            {
                const Representation& irr_a = group.getIrrep(a);
                if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;
                ntot += vrt.nalpha[a]*occ.nalpha[i];
                ntot += vrt.nbeta [a]*occ.nbeta [i];
                assert(count++ < nirrep);
            }
            assert(i < nirrep-1 || count == nirrep);
        }

        DistMatrix<U> H_elem(ntot, ntot);
        DistMatrix<U> C_elem(ntot, ntot);
        DistMatrix<U,VC,STAR> C_local;
        DistMatrix<U> E_elem;
        DistMatrix<U,STAR,STAR> E_local;

        int offbj = 0;
        for (int spin_bj : {1,0})
        {
            for (int j = 0;j < nirrep;j++)
            {
                const Representation& irr_j = group.getIrrep(j);
                for (int b = 0;b < nirrep;b++)
                {
                    const Representation& irr_b = group.getIrrep(b);
                    if (!(irr_b*irr_j*irr_R).isTotallySymmetric()) continue;

                    int nb = (spin_bj == 1 ? vrt.nalpha[b] : vrt.nbeta[b]);
                    int nj = (spin_bj == 1 ? occ.nalpha[j] : occ.nbeta[j]);
                    int nbj = nb*nj;

                    int offai = 0;
                    for (int spin_ai : {1,0})
                    {
                        for (int i = 0;i < nirrep;i++)
                        {
                            const Representation& irr_i = group.getIrrep(i);
                            for (int a = 0;a < nirrep;a++)
                            {
                                const Representation& irr_a = group.getIrrep(a);
                                if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;

                                int na = (spin_ai == 1 ? vrt.nalpha[a] : vrt.nbeta[a]);
                                int ni = (spin_ai == 1 ? occ.nalpha[i] : occ.nbeta[i]);
                                int nai = na*ni;

                                int cshift = H_elem.ColShift();
                                int rshift = H_elem.RowShift();
                                int cstride = H_elem.ColStride();
                                int rstride = H_elem.RowStride();

                                int ishift0 = cshift-offai%cstride;
                                if (ishift0 < 0) ishift0 += cstride;
                                int iloc0   = offai+ishift0;

                                int ishift1 = cshift-(offai+nai)%cstride;
                                if (ishift1 < 0) ishift1 += cstride;
                                int iloc1   = (offai+nai)+ishift1;

                                int jshift0 = rshift-offbj%rstride;
                                if (jshift0 < 0) jshift0 += rstride;
                                int jloc0   = offbj+jshift0;

                                int jshift1 = rshift-(offbj+nai)%rstride;
                                if (jshift1 < 0) jshift1 += rstride;
                                int jloc1   = (offbj+nai)+jshift1;

                                vector<tkv_pair<U>> pairs;

                                for (int iloc = iloc0;iloc < iloc1;iloc += cstride)
                                {
                                    key aidx = (iloc-offai)%na;
                                    key iidx = (iloc-offai)/na;
                                    assert(aidx >= 0 && aidx < na);
                                    assert(iidx >= 0 && iidx < ni);

                                    for (int jloc = jloc0;jloc < jloc1;jloc += rstride)
                                    {
                                        key bidx = (jloc-offbj)%nb;
                                        key jidx = (jloc-offbj)/nb;
                                        assert(bidx >= 0 && bidx < nb);
                                        assert(jidx >= 0 && jidx < nj);

                                        key k = ((iidx*nb+bidx)*nj+jidx)*na+aidx;
                                        pairs.emplace_back(k, 0);
                                    }
                                }

                                Hguess({spin_ai,spin_bj},{spin_bj,spin_ai})({a,j,b,i}).getRemoteData(pairs);

                                for (auto p : pairs)
                                {
                                    key k = p.k;
                                    int aidx = k%na;
                                    k /= na;
                                    int jidx = k%nj;
                                    k /= nj;
                                    int bidx = k%nb;
                                    k /= nb;
                                    int iidx = k;

                                    int iloc = aidx+iidx*na+offai;
                                    int jloc = bidx+jidx*nb+offbj;
                                    assert(iloc%cstride == cshift);
                                    assert(jloc%rstride == rshift);
                                    iloc /= cstride;
                                    jloc /= rstride;

                                    assert(aidx >= 0 && aidx < na);
                                    assert(bidx >= 0 && bidx < nb);
                                    assert(iidx >= 0 && iidx < ni);
                                    assert(jidx >= 0 && jidx < nj);
                                    H_elem.SetLocal(iloc, jloc, p.d);
                                }

                                offai += nai;
                            }
                        }
                    }
                    offbj += nbj;
                }
            }
        }

        HermitianEig(UPPER, H_elem, E_elem, C_elem);

        E_local = E_elem;
        C_local = C_elem;

        for (int root = 0;root < ntot;root++)
        {
            TDAevals[R].push_back(E_local.GetLocal(root, 0));
            TDAevecs[R].emplace_back("R", arena, occ.group, irr_R, vec(vrt, occ), vec(1,0), vec(0,1));
            SpinorbitalTensor<U>& evec = TDAevecs[R][root];

            vector<tkv_pair<U>> pairs;
            pairs.reserve(ntot);

            int offai = 0;
            for (int spin_ai : {1,0})
            {
                for (int i = 0;i < nirrep;i++)
                {
                    const Representation& irr_i = group.getIrrep(i);
                    for (int a = 0;a < nirrep;a++)
                    {
                        const Representation& irr_a = group.getIrrep(a);
                        if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;

                        int nai = (spin_ai == 1 ? vrt.nalpha[a] : vrt.nbeta[a])*
                                  (spin_ai == 1 ? occ.nalpha[i] : occ.nbeta[i]);

                        int cshift = C_local.ColShift();
                        int cstride = C_local.ColStride();

                        int shift0 = offai%cstride;
                        int loc0 = offai/cstride;
                        if (shift0 > cshift) loc0++;

                        int shift1 = (offai+nai)%cstride;
                        int loc1 = (offai+nai)/cstride;
                        if (shift1 > cshift) loc1++;

                        vector<tkv_pair<U>> pairs;

                        for (int loc = loc0;loc < loc1;loc++)
                        {
                            int gloc = loc*cstride+cshift-offai;
                            pairs.emplace_back(gloc, C_local.GetLocal(loc, root));
                        }

                        evec({spin_ai,0},{0,spin_ai})({a,i}).writeRemoteData(pairs);

                        offai += nai;
                    }
                }
            }
        }

        cosort(TDAevals[R].begin() , TDAevals[R].end(),
               TDAevecs[R].pbegin(), TDAevecs[R].pend());
    }

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::ElementalTDA);
REGISTER_TASK(aquarius::cc::ElementalTDA<double>, "elementaltda");
