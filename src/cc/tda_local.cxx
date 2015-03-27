#include "tda_local.hpp"

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
LocalTDA<U>::LocalTDA(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("tda.TDAevals", "TDAevals", reqs));
    addProduct(Product("tda.TDAevecs", "TDAevecs", reqs));
}

template <typename U>
bool LocalTDA<U>::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");
    const PointGroup& group = molecule.getGroup();
    int nirrep = group.getNumIrreps();

    TwoElectronOperator<U>& W = get<TwoElectronOperator<U> >("H");
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

        vector<U> data(ntot*ntot);

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

                    int nbj = (spin_bj == 1 ? vrt.nalpha[b] : vrt.nbeta[b])*
                              (spin_bj == 1 ? occ.nalpha[j] : occ.nbeta[j]);

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

                                CTFTensor<U>& this_tensor = Hguess({spin_ai,spin_bj},{spin_bj,spin_ai})({a,j,b,i});
                                CTFTensor<U> trans_tensor("trans_tensor", arena, 4, {vrt.nalpha[a],occ.nalpha[i],vrt.nbeta[b],occ.nbeta[j]}, {NS,NS,NS,NS}, true);
                                trans_tensor["ajbi"] = this_tensor["aibj"];
                                vector<U> tempdata;
                                trans_tensor.getAllData(tempdata);
                                assert(tempdata.size() == nai*nbj);
                                for (int bj = 0, aibj = 0;bj < nbj;bj++)
                                {
                                    for (int ai = 0;ai < nai;ai++, aibj++)
                                    {
                                        data[offai+ai+(offbj+bj)*ntot] = tempdata[aibj];
                                    }
                                }
                                offai += nai;
                            }
                        }
                    }
                    offbj += nbj;
                }
            }
        }

        TDAevals[R].resize(ntot);
        heev('V','U',ntot,data.data(),ntot,TDAevals[R].data());

        // cout << fixed << setprecision(5);
        // if (arena.rank == 0)
        // {
        //     cout << TDAevals[R] << endl;
        //     cout << "I'm rank 0. " << TDAevals[R][0] << endl;
        //     arena.Barrier();
        // }
        // else
        // {
        //     arena.Barrier();
        //     cout << TDAevals[R] << endl;
        //     cout << "I'm not rank 0. " << TDAevals[R][0] << endl;
        // }
        arena.Barrier();

        for (int root = 0;root < ntot;root++)
        {
            TDAevecs[R].emplace_back("R", arena, occ.group, irr_R, vec(vrt, occ), vec(1,0), vec(0,1));
            SpinorbitalTensor<U>& evec = TDAevecs[R][root];

            int offai = 0;
            for (int spin_ai = 1;spin_ai >= 0;spin_ai--)
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

                        vector<tkv_pair<U> > pairs(nai);
                        for (int ai = 0;ai < nai;ai++)
                        {
                            pairs[ai].k = ai;
                            pairs[ai].d = data[offai+ai+root*ntot];
                        }

                        if (arena.rank == 0)
                            evec({spin_ai,0},{0,spin_ai})({a,i}).writeRemoteData(pairs);
                        else
                            evec({spin_ai,0},{0,spin_ai})({a,i}).writeRemoteData();

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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LocalTDA);
REGISTER_TASK(aquarius::cc::LocalTDA<double>, "localtda");
