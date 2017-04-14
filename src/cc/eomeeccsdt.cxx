#include "util/global.hpp"

#include "convergence/davidson.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::convergence;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace cc
{

template <typename U>
class EOMEECCSDT : public Iterative<U>
{
    protected:
        Config davidson_config;
        int nroot;
        int nsinglet;
        int ntriplet;
        bool triplet = false;
        vector<U> previous;

    public:
        EOMEECCSDT(const string& name, Config& config)
        : Iterative<U>(name, config),
          davidson_config(config.get("davidson")),
          nroot(config.get<int>("nroot")),
          ntriplet(config.get<int>("ntriplet")),
          nsinglet(config.get<int>("nsinglet"))
        {
            vector<Requirement> reqs;
            reqs.emplace_back("ccsdt.T", "T");
            reqs.emplace_back("ccsdt.Hbar", "Hbar");
            reqs.emplace_back("tda.TDAevals", "TDAevals");
            reqs.emplace_back("tda.TDAevecs", "TDAevecs");
            this->addProduct("eomeeccsdt.energy", "energy", reqs);
            this->addProduct("eomeeccsdt.convergence", "convergence", reqs);
            this->addProduct("eomeeccsdt.R", "R", reqs);
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const PointGroup& group = H.getABIJ().getGroup();
            int nirrep = group.getNumIrreps();

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            this->puttmp("D", new Denominator<U>(H));
            auto& TDAevecs = this->template get<vector<unique_vector<SpinorbitalTensor<U>>>>("TDAevecs");
            auto& TDAevals = this->template get<vector<vector<U>>>("TDAevals");

            vector<tuple<U,int,int,int>> tda_sorted;

            int tot_singlet = 0;
            int tot_triplet = 0;
            for (int i = 0;i < nirrep;i++)
            {
                vector<int> spin(TDAevals[i].size());

                for (int j = 0;j < TDAevals[i].size();j++)
                {
                    auto& R1a = TDAevecs[i][j]({0,0},{0,0});
                    auto& R1b = TDAevecs[i][j]({1,0},{0,1});
                    U olap = scalar(R1a*R1b);
                    //cout << i << " " << j << " " << setprecision(15) << olap << " " << TDAevals[i][j] << endl;
                    if (olap < 0)
                    {
                        spin[j] = 1;
                        R1a -=  R1b;
                        R1b  = -R1a;
                        tot_triplet++;
                    }
                    else
                    {
                        R1a += R1b;
                        R1b  = R1a;
                        tot_singlet++;
                    }
                    U nrm = sqrt(2*scalar(R1a*R1a));
                    R1a /= nrm;
                    R1b /= nrm;
                }

                tda_sorted += zip(TDAevals[i],
                                  spin,
                                  vector<int>(TDAevals[i].size(), i),
                                  vector<int>(range<int>(TDAevals[i].size())));
            }

            assert(nroot > 0);
            assert(nsinglet >= 0);
            assert(ntriplet >= 0);
            assert(nsinglet+ntriplet <= nroot);
            assert(nsinglet <= tot_singlet);
            assert(ntriplet <= tot_triplet);
            assert(nroot <= tot_singlet+tot_triplet);

            sort(tda_sorted);

            vector<tuple<U,int,int,int>> roots;

            for (auto& root : tda_sorted)
            {
                if (nsinglet == 0) break;
                if (get<1>(root) == 0)
                {
                    nsinglet--;
                    nroot--;
                    roots.push_back(root);
                    get<1>(root) = -1;
                }
            }

            for (auto& root : tda_sorted)
            {
                if (ntriplet == 0) break;
                if (get<1>(root) == 1)
                {
                    ntriplet--;
                    nroot--;
                    roots.push_back(root);
                    get<1>(root) = -1;
                }
            }

            for (auto& root : tda_sorted)
            {
                if (nroot == 0) break;
                if (get<1>(root) != -1)
                {
                    nroot--;
                    roots.push_back(root);
                    get<1>(root) = -1;
                }
            }

            auto& Es = this->put   ("energy", new vector<U>());
            auto& Rs = this->put   ("R", new unique_vector<ExcitationOperator<U,3>>());
            auto& Cs = this->puttmp("C", new unique_vector<ExcitationOperator<U,3>>());
            auto& Zs = this->puttmp("Z", new unique_vector<ExcitationOperator<U,3>>());

            int idx = 1;
            for (int irrep = 0;irrep < nirrep;irrep++)
            {
                for (int spin : {0,1})
                {
                    triplet = spin == 1;

                    Cs.clear();
                    Zs.clear();

                    auto& davidson = this->puttmp("Davidson",
                        new Davidson<ExcitationOperator<U,3>>(davidson_config));

                    auto& rep = group.getIrrep(irrep);

                    this->puttmp(  "XMI", new SpinorbitalTensor<U>(   "X(mi)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
                    this->puttmp(  "XAE", new SpinorbitalTensor<U>(   "X(ae)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
                    this->puttmp(  "XME", new SpinorbitalTensor<U>(   "X(me)", arena, group, rep, {vrt,occ}, {0,1}, {1,0}));
                    this->puttmp("XABEJ", new SpinorbitalTensor<U>("X(ab,ej)", arena, group, rep, {vrt,occ}, {2,0}, {1,1}));
                    this->puttmp("XAMIJ", new SpinorbitalTensor<U>("X(am,ij)", arena, group, rep, {vrt,occ}, {1,1}, {0,2}));
                    this->puttmp("XABEF", new SpinorbitalTensor<U>("X(ab,ef)", arena, group, rep, {vrt,occ}, {2,0}, {2,0}));
                    this->puttmp("XAMEI", new SpinorbitalTensor<U>("X(am,ei)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
                    this->puttmp("XMNIJ", new SpinorbitalTensor<U>("X(mn,ij)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));

                    for (auto& root : roots)
                    {
                        if (spin != get<1>(root)) continue;
                        if (irrep != get<2>(root)) continue;

                        int which = get<3>(root);

                        Logger::log(arena) << "Starting root number " << idx << endl;
                        Logger::log(arena) << "Guess energy: " << fixed << setprecision(12) << get<0>(root) << endl;

                        Cs.clear();
                        Zs.clear();

                        Cs.emplace_back("R", arena, occ, vrt, group.getIrrep(irrep));
                        Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(irrep));

                        auto& R = Cs.back();

                        R(0) = 0;
                        R(1) = TDAevecs[irrep][which];
                        R(2) = 0;
                        R(3) = 0;

                        if (triplet)
                        {
                            Logger::log(arena) << "Triplet initial guess" << endl;
                        }
                        else
                        {
                            Logger::log(arena)<< "Singlet initial guess" << endl;
                        }

                        previous.assign(1, numeric_limits<U>::max());

                        Iterative<U>::run(dag, arena);

                        if (!this->isConverged())
                        {
                            this->error(arena) << "Root " << idx << " did not converge." << endl;
                        }

                        davidson.getSolution(0, R);

                        bool temp = scalar(R(1)({1,0},{0,1})*R(1)({0,0},{0,0})) < 0;
                        if (temp)
                            this->log(arena) << "triplet solution found!" << endl;
                        else
                            this->log(arena) << "singlet solution found!" << endl;
                        if (triplet != temp)
                        {
                            this->log(arena) << "WARNING: Spin character different from initial guess!" << endl;
                        }

                        R(0) = (scalar(H.getIA()["ia"]*R(1)["ai"]) +
                                0.25*scalar(H.getIJAB()["ijab"]*R(2)["abij"])) /
                               this->energy();
                        R /= sqrt(aquarius::abs(scalar(conj(R)*R)));

                        Es.push_back(this->energy());
                        Rs.emplace_back(R);

                        davidson.nextRoot();
                        idx++;
                    }
                }
            }

            this->put("convergence", new U(this->conv()));

            return true;
        }

        void iterate(const Arena& arena)
        {
            const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const SpinorbitalTensor<U>&   FME =   H.getIA();
            const SpinorbitalTensor<U>&   FAE =   H.getAB();
            const SpinorbitalTensor<U>&   FMI =   H.getIJ();
            const SpinorbitalTensor<U>& WMNEF = H.getIJAB();
            const SpinorbitalTensor<U>& WAMEF = H.getAIBC();
            const SpinorbitalTensor<U>& WABEJ = H.getABCI();
            const SpinorbitalTensor<U>& WABEF = H.getABCD();
            const SpinorbitalTensor<U>& WMNIJ = H.getIJKL();
            const SpinorbitalTensor<U>& WMNEJ = H.getIJAK();
            const SpinorbitalTensor<U>& WAMIJ = H.getAIJK();
            const SpinorbitalTensor<U>& WAMEI = H.getAIBJ();

            auto& T = this->template get<ExcitationOperator<U,3>>("T");

            auto& XMI   = this->template gettmp<SpinorbitalTensor<U>>("XMI"  );
            auto& XAE   = this->template gettmp<SpinorbitalTensor<U>>("XAE"  );
            auto& XME   = this->template gettmp<SpinorbitalTensor<U>>("XME"  );
            auto& XABEJ = this->template gettmp<SpinorbitalTensor<U>>("XABEJ");
            auto& XAMIJ = this->template gettmp<SpinorbitalTensor<U>>("XAMIJ");
            auto& XABEF = this->template gettmp<SpinorbitalTensor<U>>("XABEF");
            auto& XAMEI = this->template gettmp<SpinorbitalTensor<U>>("XAMEI");
            auto& XMNIJ = this->template gettmp<SpinorbitalTensor<U>>("XMNIJ");

            auto& D = this->template gettmp<Denominator<U>>("D");
            auto& davidson = this->template gettmp<Davidson<ExcitationOperator<U,3>>>("Davidson");

            auto& Cs = this->template gettmp<unique_vector<ExcitationOperator<U,3>>>("C");
            auto& Zs = this->template gettmp<unique_vector<ExcitationOperator<U,3>>>("Z");

            for (int root = 0;root < this->nsolution();root++)
            {
                auto& R = Cs[root];
                auto& Z = Zs[root];

                  XME[    "me"]  =      WMNEF["mnef"]*R(1)[    "fn"];

                  XMI[    "mi"]  =      WMNEJ["nmei"]*R(1)[    "en"];
                  XMI[    "mi"] +=  0.5*WMNEF["mnef"]*R(2)[  "efin"];

                  XAE[    "ae"]  =      WAMEF["amef"]*R(1)[    "fm"];
                  XAE[    "ae"] -=  0.5*WMNEF["mnef"]*R(2)[  "afmn"];

                XABEF[  "abef"]  =  0.5*WMNEF["mnef"]*R(2)[  "abmn"];
                XABEF[  "abef"] -=      WAMEF["amef"]*R(1)[    "bm"];

                XAMEI[  "amei"]  =      WMNEF["mnef"]*R(2)[  "afni"];
                XAMEI[  "amei"] +=      WAMEF["amef"]*R(1)[    "fi"];
                XAMEI[  "amei"] -=      WMNEJ["nmei"]*R(1)[    "an"];

                XMNIJ[  "mnij"]  =  0.5*WMNEF["mnef"]*R(2)[  "efij"];
                XMNIJ[  "mnij"] +=      WMNEJ["mnej"]*R(1)[    "ei"];

                XAMIJ[  "amij"]  =  0.5*WAMEF["amef"]*R(2)[  "efij"];
                XAMIJ[  "amij"] +=      WMNEJ["nmej"]*R(2)[  "aein"];
                XAMIJ[  "amij"] +=      WAMEI["amej"]*R(1)[    "ei"];
                XAMIJ[  "amij"] -=      WMNIJ["nmij"]*R(1)[    "an"];
                XAMIJ[  "amij"] +=        XME[  "me"]*T(2)[  "aeij"];
                XAMIJ[  "amij"] +=  0.5*WMNEF["mnef"]*R(3)["aefijn"];

                XABEJ[  "abej"]  =  0.5*WMNEJ["mnej"]*R(2)[  "abmn"];
                XABEJ[  "abej"] +=      WAMEF["amef"]*R(2)[  "fbmj"];
                XABEJ[  "abej"] +=      WABEF["abef"]*R(1)[    "fj"];
                XABEJ[  "abej"] -=      WAMEI["amej"]*R(1)[    "bm"];
                XABEJ[  "abej"] -=  0.5*WMNEF["mnef"]*R(3)["afbmnj"];

                 Z(1)[    "ai"]  =        FAE[  "ae"]*R(1)[    "ei"];
                 Z(1)[    "ai"] -=        FMI[  "mi"]*R(1)[    "am"];
                 Z(1)[    "ai"] -=      WAMEI["amei"]*R(1)[    "em"];
                 Z(1)[    "ai"] +=        FME[  "me"]*R(2)[  "aeim"];
                 Z(1)[    "ai"] +=  0.5*WAMEF["amef"]*R(2)[  "efim"];
                 Z(1)[    "ai"] -=  0.5*WMNEJ["mnei"]*R(2)[  "eamn"];
                 Z(1)[    "ai"] += 0.25*WMNEF["mnef"]*R(3)["efamni"];

                 Z(2)[  "abij"]  =      WABEJ["abej"]*R(1)[    "ei"];
                 Z(2)[  "abij"] -=      WAMIJ["amij"]*R(1)[    "bm"];
                 Z(2)[  "abij"] +=        FAE[  "ae"]*R(2)[  "ebij"];
                 Z(2)[  "abij"] -=        FMI[  "mi"]*R(2)[  "abmj"];
                 Z(2)[  "abij"] +=        XAE[  "ae"]*T(2)[  "ebij"];
                 Z(2)[  "abij"] -=        XMI[  "mi"]*T(2)[  "abmj"];
                 Z(2)[  "abij"] +=  0.5*WMNIJ["mnij"]*R(2)[  "abmn"];
                 Z(2)[  "abij"] +=  0.5*WABEF["abef"]*R(2)[  "efij"];
                 Z(2)[  "abij"] -=      WAMEI["amei"]*R(2)[  "ebmj"];
                 Z(2)[  "abij"] +=        FME[  "me"]*R(3)["abeijm"];
                 Z(2)[  "abij"] +=        XME[  "me"]*T(3)["abeijm"];
                 Z(2)[  "abij"] +=  0.5*WAMEF["amef"]*R(3)["efbimj"];
                 Z(2)[  "abij"] -=  0.5*WMNEJ["mnej"]*R(3)["aebimn"];

                  XMI[    "mi"] +=        FME[  "me"]*R(1)[    "ei"];

                  XAE[    "ae"] -=        FME[  "me"]*R(1)[    "am"];

                 Z(3)["abcijk"]  =      WABEJ["abej"]*R(2)[  "ecik"];
                 Z(3)["abcijk"] -=      WAMIJ["amij"]*R(2)[  "bcmk"];
                 Z(3)["abcijk"] +=      XABEJ["abej"]*T(2)[  "ecik"];
                 Z(3)["abcijk"] -=      XAMIJ["amij"]*T(2)[  "bcmk"];
                 Z(3)["abcijk"] +=        FAE[  "ae"]*R(3)["ebcijk"];
                 Z(3)["abcijk"] -=        FMI[  "mi"]*R(3)["abcmjk"];
                 Z(3)["abcijk"] +=        XAE[  "ae"]*T(3)["ebcijk"];
                 Z(3)["abcijk"] -=        XMI[  "mi"]*T(3)["abcmjk"];
                 Z(3)["abcijk"] +=  0.5*WABEF["abef"]*R(3)["efcijk"];
                 Z(3)["abcijk"] -=      WAMEI["amei"]*R(3)["ebcmjk"];
                 Z(3)["abcijk"] +=  0.5*WMNIJ["mnij"]*R(3)["abcmnk"];
                 Z(3)["abcijk"] +=  0.5*XABEF["abef"]*T(3)["efcijk"];
                 Z(3)["abcijk"] -=      XAMEI["amei"]*T(3)["ebcmjk"];
                 Z(3)["abcijk"] +=  0.5*XMNIJ["mnij"]*T(3)["abcmnk"];
            }

            vector<U> energies = davidson.extrapolate(Cs, Zs, D);

            for (int i = 0;i < this->nsolution();i++)
            {
                this->energy(i) = energies[i];
                this->conv(i) = Zs[i].norm(00);
            }
        }
};

}
}

static const char* spec = R"(

nroot?
    int 1,
nsinglet?
    int 0,
ntriplet?
    int 0,
convergence?
    double 1e-9,
max_iterations?
    int 150,
conv_type?
    enum { MAXE, RMSE, MAE },
davidson?
{
    damping?
            double 0.0,
    start?
            int 1,
    order?
            int 10,
    jacobi?
            bool false,
    num_reduce?
            int 3,
    compaction?
            enum { discrete, continuous },
}

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSDT);
REGISTER_TASK(aquarius::cc::EOMEECCSDT<double>, "eomeeccsdt",spec);
