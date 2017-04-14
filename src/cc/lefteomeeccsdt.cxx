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
class LeftEOMEECCSDT : public Iterative<U>
{
    protected:
        Config davidson_config;
        bool triplet = false;
        vector<U> previous;

    public:
        LeftEOMEECCSDT(const string& name, Config& config)
        : Iterative<U>(name, config),
          davidson_config(config.get("davidson"))
        {
            vector<Requirement> reqs;
            reqs.emplace_back("ccsdt.T", "T");
            reqs.emplace_back("ccsdt.Hbar", "Hbar");
            reqs.emplace_back("eomeeccsdt.energy", "energy");
            reqs.emplace_back("eomeeccsdt.R", "R");
            this->addProduct("eomeeccsdt.L", "L", reqs);
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const PointGroup& group = H.getABIJ().getGroup();
            int nirrep = group.getNumIrreps();

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            this->puttmp("D", new Denominator<U>(H));

            auto& Es = this->template get<vector<U>>("energy");
            auto& Rs = this->template get<unique_vector<ExcitationOperator<U,3>>>("R");

            auto& Ls = this->put   ("L", new unique_vector<DeexcitationOperator<U,3>>());
            auto& Cs = this->puttmp("C", new unique_vector<DeexcitationOperator<U,3>>());
            auto& Zs = this->puttmp("Z", new unique_vector<DeexcitationOperator<U,3>>());

            int idx = 1;
            for (int irrep = 0;irrep < nirrep;irrep++)
            {
                Cs.clear();
                Zs.clear();

                auto& davidson = this->puttmp("Davidson",
                    new Davidson<DeexcitationOperator<U,3>>(davidson_config));

                auto& rep = group.getIrrep(irrep);

                this->puttmp(  "XEM", new SpinorbitalTensor<U>(   "X(em)", arena, group, rep, {vrt,occ}, {1,0}, {0,1}));
                this->puttmp(  "XIM", new SpinorbitalTensor<U>(   "X(im)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
                this->puttmp(  "XEA", new SpinorbitalTensor<U>(   "X(ea)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
                this->puttmp("XEJAB", new SpinorbitalTensor<U>("X(ej,ab)", arena, group, rep, {vrt,occ}, {1,1}, {2,0}));
                this->puttmp("XIJAM", new SpinorbitalTensor<U>("X(ij,am)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));
                this->puttmp("XEFAB", new SpinorbitalTensor<U>("X(ef,ab)", arena, group, rep, {vrt,occ}, {2,0}, {2,0}));
                this->puttmp("XEIAM", new SpinorbitalTensor<U>("X(ei,am)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
                this->puttmp("XIJMN", new SpinorbitalTensor<U>("X(ij,mn)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));

                for (int root = 0;root < Rs.size();root++)
                {
                    auto& R = Rs[root];
                    U E = Es[root];

                    if (!(R.getRepresentation()*group.getIrrep(irrep))) continue;

                    davidson.reset(E);

                    Logger::log(arena) << "Starting root number " << idx << endl;
                    Logger::log(arena) << "Guess energy: " << fixed << setprecision(12) << E << endl;

                    Cs.clear();
                    Zs.clear();

                    Cs.emplace_back("C", arena, occ, vrt, group.getIrrep(irrep));
                    Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(irrep));

                    auto& L = Cs.back();
                    L(0)           =              0;
                    L(1)[    "ia"] = R(1)[    "ai"];
                    L(2)[  "ijab"] = R(2)[  "abij"];
                    L(3)["ijkabc"] = R(3)["abcijk"];
                    L /= sqrt(aquarius::abs(scalar(conj(L)*L)));

                    triplet = scalar(L(1)({0,1},{1,0})*L(1)({0,0},{0,0})) < 0;

                    previous.assign(1, numeric_limits<U>::max());

                    Iterative<U>::run(dag, arena);

                    if (!this->isConverged())
                    {
                        this->error(arena) << "Root " << idx << " did not converge." << endl;
                    }

                    davidson.getSolution(0, L);

                    bool temp = scalar(L(1)({0,1},{1,0})*L(1)({0,0},{0,0})) < 0;
                    if (temp)
                        this->log(arena) << "triplet solution found!" << endl;
                    else
                        this->log(arena) << "singlet solution found!" << endl;
                    if (triplet != temp)
                    {
                        this->log(arena) << "WARNING: Spin character different from initial guess!" << endl;
                    }

                    L /=  (1.0/ 1.0)*scalar(conj(L(1)[    "ia"])*R(1)[    "ai"]) +
                          (1.0/ 4.0)*scalar(conj(L(2)[  "ijab"])*R(2)[  "abij"]) +
                          (1.0/36.0)*scalar(conj(L(3)["ijkabc"])*R(3)["abcijk"]);

                    Ls.emplace_back(L);

                    davidson.nextRoot();
                    idx++;
                }
            }

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

            auto& XEM   = this->template gettmp<SpinorbitalTensor<U>>(  "XEM");
            auto& XIM   = this->template gettmp<SpinorbitalTensor<U>>(  "XIM");
            auto& XEA   = this->template gettmp<SpinorbitalTensor<U>>(  "XEA");
            auto& XEJAB = this->template gettmp<SpinorbitalTensor<U>>("XEJAB");
            auto& XIJAM = this->template gettmp<SpinorbitalTensor<U>>("XIJAM");
            auto& XEFAB = this->template gettmp<SpinorbitalTensor<U>>("XEFAB");
            auto& XEIAM = this->template gettmp<SpinorbitalTensor<U>>("XEIAM");
            auto& XIJMN = this->template gettmp<SpinorbitalTensor<U>>("XIJMN");

            auto& D = this->template gettmp<Denominator<U>>("D");
            auto& davidson = this->template gettmp<Davidson<DeexcitationOperator<U,3>>>("Davidson");

            auto& Cs = this->template gettmp<unique_vector<DeexcitationOperator<U,3>>>("C");
            auto& Zs = this->template gettmp<unique_vector<DeexcitationOperator<U,3>>>("Z");

            for (int root = 0;root < this->nsolution();root++)
            {
                auto& L = Cs[root];
                auto& Z = Zs[root];

                Z(1) = 0;

                  XIM[    "im"]  =  (1.0/12.0)*T(3)["efgmno"]* L(3)["inoefg"];
                  XEA[    "ea"]  = -(1.0/12.0)*T(3)["efgmno"]* L(3)["mnoafg"];

                 Z(1)[    "ia"]  =              FME[    "ie"]*  XEA[    "ea"];
                 Z(1)[    "ia"] -=              FME[    "ma"]*  XIM[    "im"];

                XEFAB[  "efab"]  =  (1.0/ 6.0)*T(3)["efgmno"]* L(3)["mnoabg"];
                XEIAM[  "eiam"]  = -(1.0/ 4.0)*T(3)["efgmno"]* L(3)["inoafg"];
                XIJMN[  "ijmn"]  =  (1.0/ 6.0)*T(3)["efgmno"]* L(3)["ijoefg"];

                XIJAM[  "ijam"]  =         0.5*T(2)[  "efmn"]* L(3)["ijnaef"];
                XEJAB[  "ejab"]  = -       0.5*T(2)[  "efmn"]* L(3)["mjnabf"];

                  XEM[    "em"]  =        0.25*T(3)["efgmno"]* L(2)[  "nofg"];
                  XEM[    "em"] -=         0.5*T(2)[  "efno"]*XIJAM[  "onfm"];

                  XIM[    "ij"] +=         0.5*T(2)[  "efjm"]* L(2)[  "imef"];
                  XEA[    "ab"] -=         0.5*T(2)[  "aemn"]* L(2)[  "mnbe"];

                 Z(1)[    "ia"] +=              FAE[    "ea"]* L(1)[    "ie"];
                 Z(1)[    "ia"] -=              FMI[    "im"]* L(1)[    "ma"];
                 Z(1)[    "ia"] -=            WAMEI[  "eiam"]* L(1)[    "me"];
                 Z(1)[    "ia"] +=        0.5*WABEJ[  "efam"]* L(2)[  "imef"];
                 Z(1)[    "ia"] -=        0.5*WAMIJ[  "eimn"]* L(2)[  "mnea"];
                 Z(1)[    "ia"] -=            WMNEJ[  "inam"]*  XIM[    "mn"];
                 Z(1)[    "ia"] -=            WAMEF[  "fiea"]*  XEA[    "ef"];
                 Z(1)[    "ia"] +=            WMNEF[  "miea"]*  XEM[    "em"];
                 Z(1)[    "ia"] -=        0.5*WABEF[  "efga"]*XEJAB[  "gief"];
                 Z(1)[    "ia"] +=            WAMEI[  "eifm"]*XEJAB[  "fmea"];
                 Z(1)[    "ia"] -=            WAMEI[  "eman"]*XIJAM[  "inem"];
                 Z(1)[    "ia"] +=        0.5*WMNIJ[  "imno"]*XIJAM[  "noam"];
                 Z(1)[    "ia"] -=        0.5*WAMEF[  "gief"]*XEFAB[  "efga"];
                 Z(1)[    "ia"] +=            WAMEF[  "fmea"]*XEIAM[  "eifm"];
                 Z(1)[    "ia"] -=            WMNEJ[  "inem"]*XEIAM[  "eman"];
                 Z(1)[    "ia"] +=        0.5*WMNEJ[  "noam"]*XIJMN[  "imno"];

                 Z(2)[  "ijab"]  =              FME[    "ia"]* L(1)[    "jb"];
                 Z(2)[  "ijab"] +=            WAMEF[  "ejab"]* L(1)[    "ie"];
                 Z(2)[  "ijab"] -=            WMNEJ[  "ijam"]* L(1)[    "mb"];
                 Z(2)[  "ijab"] +=              FAE[    "ea"]* L(2)[  "ijeb"];
                 Z(2)[  "ijab"] -=              FMI[    "im"]* L(2)[  "mjab"];
                 Z(2)[  "ijab"] +=        0.5*WABEF[  "efab"]* L(2)[  "ijef"];
                 Z(2)[  "ijab"] +=        0.5*WMNIJ[  "ijmn"]* L(2)[  "mnab"];
                 Z(2)[  "ijab"] +=            WAMEI[  "eiam"]* L(2)[  "mjbe"];
                 Z(2)[  "ijab"] -=            WMNEF[  "mjab"]*  XIM[    "im"];
                 Z(2)[  "ijab"] +=            WMNEF[  "ijeb"]*  XEA[    "ea"];
                 Z(2)[  "ijab"] +=        0.5*WMNEF[  "ijef"]*XEFAB[  "efab"];
                 Z(2)[  "ijab"] +=            WMNEF[  "imea"]*XEIAM[  "ejbm"];
                 Z(2)[  "ijab"] +=        0.5*WMNEF[  "mnab"]*XIJMN[  "ijmn"];
                 Z(2)[  "ijab"] -=            WAMEF[  "fiae"]*XEJAB[  "ejbf"];
                 Z(2)[  "ijab"] -=            WMNEJ[  "ijem"]*XEJAB[  "emab"];
                 Z(2)[  "ijab"] -=            WAMEF[  "emab"]*XIJAM[  "ijem"];
                 Z(2)[  "ijab"] -=            WMNEJ[  "niam"]*XIJAM[  "mjbn"];
                 Z(2)[  "ijab"] +=        0.5*WABEJ[  "efbm"]* L(3)["ijmaef"];
                 Z(2)[  "ijab"] -=        0.5*WAMIJ[  "ejnm"]* L(3)["imnabe"];

                 Z(3)["ijkabc"]  =            WMNEF[  "ijab"]* L(1)[    "kc"];
                 Z(3)["ijkabc"] +=              FME[    "ia"]* L(2)[  "jkbc"];
                 Z(3)["ijkabc"] +=            WAMEF[  "ekbc"]* L(2)[  "ijae"];
                 Z(3)["ijkabc"] -=            WMNEJ[  "ijam"]* L(2)[  "mkbc"];
                 Z(3)["ijkabc"] +=            WMNEF[  "ijae"]*XEJAB[  "ekbc"];
                 Z(3)["ijkabc"] -=            WMNEF[  "mkbc"]*XIJAM[  "ijam"];
                 Z(3)["ijkabc"] +=              FAE[    "ea"]* L(3)["ijkebc"];
                 Z(3)["ijkabc"] -=              FMI[    "im"]* L(3)["mjkabc"];
                 Z(3)["ijkabc"] +=        0.5*WABEF[  "efab"]* L(3)["ijkefc"];
                 Z(3)["ijkabc"] +=        0.5*WMNIJ[  "ijmn"]* L(3)["mnkabc"];
                 Z(3)["ijkabc"] +=            WAMEI[  "eiam"]* L(3)["mjkbec"];

                 //printf("Z1: %.15f\n", Z(1)({0,1},{1,0}).norm(2));
                 //printf("Z2: %.15f\n", Z(2)({0,1},{1,0}).norm(2));
                 //printf("Z3: %.15f\n", (1.0/36.0)*scalar(Z(3)*Z(3)));
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LeftEOMEECCSDT);
REGISTER_TASK(aquarius::cc::LeftEOMEECCSDT<double>, "lefteomeeccsdt",spec);
