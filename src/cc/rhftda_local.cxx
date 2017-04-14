#include "util/global.hpp"

#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"
#include "input/molecule.hpp"

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
class LocalRHFTDA : public Task
{
    public:
        LocalRHFTDA(const string& name, Config& config)
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

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const Molecule& molecule = this->template get<Molecule>("molecule");
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

                vector<U> data_sing(ntot*ntot);
                vector<U> data_trip(ntot*ntot);

                int offbj = 0;
                for (int j = 0;j < nirrep;j++)
                {
                    const Representation& irr_j = group.getIrrep(j);
                    for (int b = 0;b < nirrep;b++)
                    {
                        const Representation& irr_b = group.getIrrep(b);
                        if (!(irr_b*irr_j*irr_R).isTotallySymmetric()) continue;

                        int nbj = nA[b]*nI[j];

                        int offai = 0;
                        for (int i = 0;i < nirrep;i++)
                        {
                            const Representation& irr_i = group.getIrrep(i);
                            for (int a = 0;a < nirrep;a++)
                            {
                                const Representation& irr_a = group.getIrrep(a);
                                if (!(irr_a*irr_i*irr_R).isTotallySymmetric()) continue;

                                int nai = nA[a]*nI[i];

                                vector<U> data;

                                H_sing({a,i,b,j}).getAllData(data);
                                assert(data.size() == nai*nbj);
                                for (int bj = 0;bj < nbj;bj++)
                                {
                                    for (int ai = 0;ai < nai;ai++)
                                    {
                                        data_sing[offai+ai+(offbj+bj)*ntot] = data[ai+bj*nai];
                                    }
                                }

                                H_trip({a,i,b,j}).getAllData(data);
                                assert(data.size() == nai*nbj);
                                for (int bj = 0;bj < nbj;bj++)
                                {
                                    for (int ai = 0;ai < nai;ai++)
                                    {
                                        data_trip[offai+ai+(offbj+bj)*ntot] = data[ai+bj*nai];
                                    }
                                }

                                offai += nai;
                            }
                        }
                        offbj += nbj;
                    }
                }

                evals_sing[R].resize(ntot);
                evals_trip[R].resize(ntot);
                heev('V','U',ntot,data_sing.data(),ntot,evals_sing[R].data());
                heev('V','U',ntot,data_trip.data(),ntot,evals_trip[R].data());

                arena.comm().Barrier();

                for (int root = 0;root < ntot;root++)
                {
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

                            vector<tkv_pair<U>> pairs(nai);

                            for (int ai = 0;ai < nai;ai++)
                            {
                                pairs[ai].k = ai;
                                pairs[ai].d = data_sing[offai+ai+root*ntot];
                            }

                            if (arena.rank == 0) evec_sing({a,i}).writeRemoteData(pairs);
                            else                 evec_sing({a,i}).writeRemoteData();

                            for (int ai = 0;ai < nai;ai++)
                            {
                                pairs[ai].k = ai;
                                pairs[ai].d = data_trip[offai+ai+root*ntot];
                            }

                            if (arena.rank == 0) evec_trip({a,i}).writeRemoteData(pairs);
                            else                 evec_trip({a,i}).writeRemoteData();

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
};

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LocalRHFTDA);
REGISTER_TASK(aquarius::cc::LocalRHFTDA<double>, "localrhftda");
