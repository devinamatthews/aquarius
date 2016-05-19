#include "rhftda_elemental.hpp"

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
ElementalRHFTDA<U>::ElementalRHFTDA(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    reqs.push_back(Requirement("mofock", "f"));
    reqs.push_back(Requirement("<Ai|Bj>", "VAIBJ"));
    reqs.push_back(Requirement("<Ai|Jb>", "VAIJB"));
    addProduct(Product("rhftda.evals_sing", "evals_sing", reqs));
    addProduct(Product("rhftda.evals_trip", "evals_trip", reqs));
    addProduct(Product("rhftda.evecs_sing", "evecs_sing", reqs));
    addProduct(Product("rhftda.evecs_trip", "evecs_trip", reqs));
}

template <typename U>
bool ElementalRHFTDA<U>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    int nirrep = group.getNumIrreps();

    const auto& f = get<OneElectronOperator<U>>("f");
    const Space& occ = f.occ;
    const Space& vrt = f.vrt;
    const vector<int>& nI = occ.nalpha;
    const vector<int>& nA = vrt.nalpha;

    const auto& VAIBJ = this->template get<SymmetryBlockedTensor<U>>("VAIBJ");
    const auto& VAIJB = this->template get<SymmetryBlockedTensor<U>>("VAIJB");
    const auto& fAB = f.getAB()({0,0},{0,0});
    const auto& fIJ = f.getIJ()({0,0},{0,0});

    SymmetryBlockedTensor<U> H_sing("H_sing", arena, group, 4, {nA,nI,nA,nI}, {NS,NS,NS,NS}, true);
    SymmetryBlockedTensor<U> H_trip("H_trip", arena, group, 4, {nA,nI,nA,nI}, {NS,NS,NS,NS}, true);

    H_trip["aibi"]  =     fAB[  "ab"];
    H_trip["aibj"] -=   VAIBJ["ajbi"];
    H_trip["aiaj"] -=     fIJ[  "ji"];

    H_sing["aibj"]  =  H_trip["aibj"];
    H_sing["aibj"] += 2*VAIJB["ajib"];

    auto& evecs_sing = put("evecs_sing", new vector<unique_vector<SymmetryBlockedTensor<U>>>(nirrep));
    auto& evecs_trip = put("evecs_trip", new vector<unique_vector<SymmetryBlockedTensor<U>>>(nirrep));
    auto& evals_sing = put("evals_sing", new vector<vector<U>>(nirrep));
    auto& evals_trip = put("evals_trip", new vector<vector<U>>(nirrep));

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
                ntot += nA[a]*nI[i];
                assert(count++ < nirrep);
            }
            assert(i < nirrep-1 || count == nirrep);
        }

        DistMatrix<U> H_sing_elem(ntot, ntot);
        DistMatrix<U> C_sing_elem(ntot, ntot);
        DistMatrix<U> H_trip_elem(ntot, ntot);
        DistMatrix<U> C_trip_elem(ntot, ntot);
        DistMatrix<U,VC,STAR> C_sing_local;
        DistMatrix<U,VC,STAR> C_trip_local;
        DistMatrix<U> E_sing_elem;
        DistMatrix<U> E_trip_elem;
        DistMatrix<U,STAR,STAR> E_sing_local;
        DistMatrix<U,STAR,STAR> E_trip_local;

        int offbj = 0;
        for (int j = 0;j < nirrep;j++)
        {
            const Representation& irr_j = group.getIrrep(j);
            for (int b = 0;b < nirrep;b++)
            {
                const Representation& irr_b = group.getIrrep(b);
                if (!(irr_b*irr_j*irr_R).isTotallySymmetric()) continue;

                int nb = nA[b];
                int nj = nI[j];
                int nbj = nb*nj;

                int offai = 0;
                for (int i = 0;i < nirrep;i++)
                {
                    const Representation& irr_i = group.getIrrep(i);
                    for (int a = 0;a < nirrep;a++)
                    {
                        const Representation& irr_a = group.getIrrep(a);
                        if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;

                        int na = nA[a];
                        int ni = nI[i];
                        int nai = na*ni;

                        int cshift = H_sing_elem.ColShift();
                        int rshift = H_sing_elem.RowShift();
                        int cstride = H_sing_elem.ColStride();
                        int rstride = H_sing_elem.RowStride();

                        int ishift0 = cshift-offai%cstride;
                        if (ishift0 < 0) ishift0 += cstride;
                        int iloc0   = offai+ishift0;
                        assert(iloc0%cstride == cshift);

                        int ishift1 = cshift-(offai+nai)%cstride;
                        if (ishift1 < 0) ishift1 += cstride;
                        int iloc1   = (offai+nai)+ishift1;
                        assert(iloc1%cstride == cshift);

                        int jshift0 = rshift-offbj%rstride;
                        if (jshift0 < 0) jshift0 += rstride;
                        int jloc0   = offbj+jshift0;
                        assert(jloc0%rstride == rshift);

                        int jshift1 = rshift-(offbj+nai)%rstride;
                        if (jshift1 < 0) jshift1 += rstride;
                        int jloc1   = (offbj+nai)+jshift1;
                        assert(jloc1%rstride == rshift);

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

                                key k = ((jidx*nb+bidx)*ni+iidx)*na+aidx;
                                pairs.emplace_back(k, 0);
                            }
                        }

                        H_sing({a,i,b,j}).getRemoteData(pairs);

                        for (auto p : pairs)
                        {
                            key k = p.k;
                            int aidx = k%na;
                            k /= na;
                            int iidx = k%ni;
                            k /= ni;
                            int bidx = k%nb;
                            k /= nb;
                            int jidx = k;

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
                            H_sing_elem.SetLocal(iloc, jloc, p.d);
                        }

                        H_trip({a,i,b,j}).getRemoteData(pairs);

                        for (auto p : pairs)
                        {
                            key k = p.k;
                            int aidx = k%na;
                            k /= na;
                            int iidx = k%ni;
                            k /= ni;
                            int bidx = k%nb;
                            k /= nb;
                            int jidx = k;

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
                            H_trip_elem.SetLocal(iloc, jloc, p.d);
                        }

                        offai += nai;
                    }
                }
                offbj += nbj;
            }
        }

        HermitianEig(UPPER, H_sing_elem, E_sing_elem, C_sing_elem);
        HermitianEig(UPPER, H_trip_elem, E_trip_elem, C_trip_elem);

        E_sing_local = E_sing_elem;
        E_trip_local = E_trip_elem;
        C_sing_local = C_sing_elem;
        C_trip_local = C_trip_elem;

        for (int root = 0;root < ntot;root++)
        {
            evals_sing[R].push_back(E_sing_local.GetLocal(root, 0));
            evals_trip[R].push_back(E_trip_local.GetLocal(root, 0));
            evecs_sing[R].emplace_back("R", arena, occ.group, 2, vec(nA,nI), vec((int)NS,(int)NS), true);
            evecs_trip[R].emplace_back("R", arena, occ.group, 2, vec(nA,nI), vec((int)NS,(int)NS), true);
            auto& evec_sing = evecs_sing[R][root];
            auto& evec_trip = evecs_trip[R][root];

            int offai = 0;
            for (int i = 0;i < nirrep;i++)
            {
                const Representation& irr_i = group.getIrrep(i);
                for (int a = 0;a < nirrep;a++)
                {
                    const Representation& irr_a = group.getIrrep(a);
                    if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;

                    int nai = nA[a]*nI[i];

                    int cshift = C_sing_local.ColShift();
                    int cstride = C_sing_local.ColStride();

                    int shift0 = offai%cstride;
                    int loc0 = offai/cstride;
                    if (shift0 > cshift) loc0++;

                    int shift1 = (offai+nai)%cstride;
                    int loc1 = (offai+nai)/cstride;
                    if (shift1 > cshift) loc1++;

                    vector<tkv_pair<U>> pairs(loc1-loc0);

                    for (int loc = loc0;loc < loc1;loc++)
                    {
                        pairs[loc-loc0].k = loc*cstride+cshift-offai;
                        pairs[loc-loc0].d = C_sing_local.GetLocal(loc, root);
                    }

                    evec_sing({a,i}).writeRemoteData(pairs);

                    for (int loc = loc0;loc < loc1;loc++)
                    {
                        pairs[loc-loc0].k = loc*cstride+cshift-offai;
                        pairs[loc-loc0].d = C_trip_local.GetLocal(loc, root);
                    }

                    evec_trip({a,i}).writeRemoteData(pairs);

                    offai += nai;
                }
            }
        }

        cosort(evals_sing[R].begin() , evals_sing[R].end(),
               evecs_sing[R].pbegin(), evecs_sing[R].pend());
        cosort(evals_trip[R].begin() , evals_trip[R].end(),
               evecs_trip[R].pbegin(), evecs_trip[R].pend());
    }

    return true;
}

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::ElementalRHFTDA);
REGISTER_TASK(aquarius::cc::ElementalRHFTDA<double>, "elementalrhftda");
