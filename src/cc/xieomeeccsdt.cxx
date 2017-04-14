#include "util/global.hpp"

#include "operator/st2eoperator.hpp"
#include "operator/deexcitationoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "convergence/diis.hpp"
#include "util/iterative.hpp"
#include "task/task.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::convergence;

namespace aquarius
{
namespace cc
{

template <typename U>
class XiEOMEECCSDT : public Task
{
    public:
        XiEOMEECCSDT(const string& name, Config& config)
        : Task(name, config)
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("ccsdt.Hbar", "Hbar"));
            reqs.push_back(Requirement("eomeeccsdt.R", "R"));
            reqs.push_back(Requirement("eomeeccsdt.L", "L"));
            this->addProduct(Product("eomeeccsdt.Xi", "Xi", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;
            auto& group = occ.group;

            auto& Xis = this->put("Xi", new unique_vector<DeexcitationOperator<U,3>>());

            auto& Rs = this->template get<unique_vector<ExcitationOperator<U,3>>>("R");
            auto& Ls = this->template get<unique_vector<ExcitationOperator<U,3>>>("L");

            auto&   FME =   H.getIA();
            auto&   FAE =   H.getAB();
            auto&   FMI =   H.getIJ();
            auto& WMNEF = H.getIJAB();
            auto& WAMEF = H.getAIBC();
            auto& WABEJ = H.getABCI();
            auto& WABEF = H.getABCD();
            auto& WMNIJ = H.getIJKL();
            auto& WMNEJ = H.getIJAK();
            auto& WAMIJ = H.getAIJK();
            auto& WAMEI = H.getAIBJ();

            for (int root = 0;root < Rs.size();root++)
            {
                Xis.emplace_back("Xi", arena, occ, vrt);
                auto& Xi = Xis.back();

                auto& R = Rs[root];
                auto& L = Ls[root];
                auto& rep = R.getRepresentation();

                auto&   DIA = this->puttmp(  "DIA", new SpinorbitalTensor<U>(   "D(ia)", arena, group, rep, {vrt,occ}, {0,1}, {1,0}));
                auto&   DAI = this->puttmp(  "DAI", new SpinorbitalTensor<U>(   "D(ai)", arena, group, rep, {vrt,occ}, {1,0}, {0,1}));
                auto&   DAB = this->puttmp(  "DAB", new SpinorbitalTensor<U>(   "D(ab)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
                auto&   DIJ = this->puttmp(  "DIJ", new SpinorbitalTensor<U>(   "D(ij)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
                auto& GIJKL = this->puttmp("GIJKL", new SpinorbitalTensor<U>("G(ij,kl)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));
                auto& GAIBJ = this->puttmp("GAIBJ", new SpinorbitalTensor<U>("G(ai,bj)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
                auto& GABCD = this->puttmp("GABCD", new SpinorbitalTensor<U>("G(ab,cd)", arena, group, rep, {vrt,occ}, {2,0}, {2,0}));
                auto& GAIBC = this->puttmp("GAIBC", new SpinorbitalTensor<U>("G(ai,bc)", arena, group, rep, {vrt,occ}, {1,1}, {2,0}));
                auto& GIJAK = this->puttmp("GIJAK", new SpinorbitalTensor<U>("G(ij,ak)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));
                auto& GIJAB = this->puttmp("GIJAB", new SpinorbitalTensor<U>("G(ij,ab)", arena, group, rep, {vrt,occ}, {0,2}, {2,0}));

                auto&   XIA = this->puttmp(  "XIA", new SpinorbitalTensor<U>(   "X(ia)", arena, group, rep, {vrt,occ}, {0,1}, {1,0}));
                auto&   XEM = this->puttmp(  "XEM", new SpinorbitalTensor<U>(   "X(em)", arena, group, rep, {vrt,occ}, {1,0}, {0,1}));
                auto&   XEA = this->puttmp(  "XEA", new SpinorbitalTensor<U>(   "X(ea)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
                auto&   XIM = this->puttmp(  "XIM", new SpinorbitalTensor<U>(   "X(im)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
                auto& XEFAB = this->puttmp("XEFAB", new SpinorbitalTensor<U>("X(ef,ab)", arena, group, rep, {vrt,occ}, {2,0}, {2,0}));
                auto& XEIAM = this->puttmp("XEIAM", new SpinorbitalTensor<U>("X(ei,am)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
                auto& XIJMN = this->puttmp("XIJMN", new SpinorbitalTensor<U>("X(ij,mn)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));
                auto& XEIMN = this->puttmp("XEIMN", new SpinorbitalTensor<U>("X(ei,mn)", arena, group, rep, {vrt,occ}, {1,1}, {0,2}));
                auto& XEFAM = this->puttmp("XEFAM", new SpinorbitalTensor<U>("X(ef,am)", arena, group, rep, {vrt,occ}, {2,0}, {1,1}));
                auto& XABIJ = this->puttmp("XABIJ", new SpinorbitalTensor<U>("X(ab,ij)", arena, group, rep, {vrt,occ}, {2,0}, {0,2}));

                U LR1 = scalar(L(1)["ia"]*R(1)["ai"]);
                U LR2 = 0.25*scalar(L(2)["ijab"]*R(2)["abij"]);
                U LR3 = U(1) - LR1 - LR2;

                Xi(0)         =                       0;
                Xi(1)[  "ia"] =      LR3 *  FME[  "ia"];
                Xi(2)["ijab"] = (LR2+LR3)*WMNEF["ijab"];
                Xi(3)         =                       0;

                if (rep)
                {
                    U FR =      scalar(  FME[  "ia"]*R(1)[  "ai"]);
                    U WR = 0.25*scalar(WMNEF["ijab"]*R(2)["abij"]);

                    Xi(2)[  "ijab"] +=     WR *L(2)[  "ijab"];
                    Xi(3)["ijkabc"] += (FR+WR)*L(3)["ijkabc"];
                }

                  XEM[    "em"]  =  (1.0/ 4.0)*WMNEF[  "nofg"]* R(3)["efgmno"];

                  DAI[    "ai"]  =  (1.0/ 4.0)* L(2)[  "mnef"]* R(3)["aefimn"];

                  DIJ[    "ij"]  =  (1.0/12.0)* L(3)["imnefg"]* R(3)["efgjmn"];

                  DAB[    "ab"]  = -(1.0/12.0)* L(3)["mnobef"]* R(3)["aefmno"];

                GABCD[  "efab"]  =  (1.0/ 6.0)* L(3)["mnoabg"]* R(3)["efgmno"];

                GAIBJ[  "eiam"]  = -(1.0/ 4.0)* L(3)["inoafg"]* R(3)["efgmno"];

                GIJKL[  "ijmn"]  =  (1.0/ 6.0)* L(3)["ijoefg"]* R(3)["efgmno"];

                XEIMN[  "eimn"]  =  (1.0/ 2.0)*WMNEF[  "iofg"]* R(3)["efgmno"];

                XEFAM[  "efam"]  = -(1.0/ 2.0)*WMNEF[  "noag"]* R(3)["fegmno"];

                XABIJ[  "abij"]  =               FME[    "me"]* R(3)["abeijm"];
                XABIJ[  "abij"] +=  (1.0/ 2.0)*WAMEF[  "amef"]* R(3)["ebfijm"];
                XABIJ[  "abij"] -=  (1.0/ 2.0)*WMNEJ[  "mnej"]* R(3)["abeinm"];

                Xi(1)[    "ia"] +=               FME[    "ie"]*  DAB[    "ea"];
                Xi(1)[    "ia"] -=               FME[    "ma"]*  DIJ[    "im"];
                Xi(1)[    "ia"] +=             WMNEF[  "imae"]*  DAI[    "em"];
                Xi(1)[    "ia"] -=             WMNEJ[  "iman"]*  DIJ[    "nm"];
                Xi(1)[    "ia"] -=             WAMEF[  "eifa"]*  DAB[    "fe"];
                Xi(1)[    "ia"] +=               XEM[    "em"]* L(2)[  "imae"];
                Xi(1)[    "ia"] +=  (1.0/ 2.0)*XEFAM[  "efam"]* L(2)[  "imef"];
                Xi(1)[    "ia"] -=  (1.0/ 2.0)*XEIMN[  "eimn"]* L(2)[  "mnea"];
                Xi(1)[    "ia"] -=  (1.0/ 2.0)*WAMEF[  "gief"]*GABCD[  "efga"];
                Xi(1)[    "ia"] +=             WAMEF[  "fmea"]*GAIBJ[  "eifm"];
                Xi(1)[    "ia"] -=             WMNEJ[  "inem"]*GAIBJ[  "eman"];
                Xi(1)[    "ia"] +=  (1.0/ 2.0)*WMNEJ[  "noam"]*GIJKL[  "imno"];
                Xi(1)[    "ia"] +=  (1.0/ 4.0)*XABIJ[  "efmn"]* L(3)["imnaef"];

                  XEM[    "em"] +=               FME[    "nf"]* R(2)[  "efmn"];
                  XEM[    "em"] +=  (1.0/ 2.0)*WAMEF[  "enfg"]* R(2)[  "fgmn"];
                  XEM[    "em"] -=  (1.0/ 2.0)*WMNEJ[  "onfm"]* R(2)[  "efno"];

                  XEA[    "ea"]  = -(1.0/ 2.0)*WMNEF[  "mnaf"]* R(2)[  "efmn"];

                  XIM[    "im"]  =  (1.0/ 2.0)*WMNEF[  "inef"]* R(2)[  "efmn"];

                  DIA[    "ia"]  =  (1.0/ 4.0)* L(3)["imnaef"]* R(2)[  "efmn"];

                  DIJ[    "ij"] +=  (1.0/ 2.0)* L(2)[  "imef"]* R(2)[  "efjm"];

                  DAB[    "ab"] -=  (1.0/ 2.0)* L(2)[  "mnbe"]* R(2)[  "aemn"];

                GAIBC[  "aibc"]  = -(1.0/ 2.0)* L(3)["minbce"]* R(2)[  "aemn"];

                GIJAK[  "ijak"]  =  (1.0/ 2.0)* L(3)["ijmaef"]* R(2)[  "efkm"];

                XEFAB[  "efab"]  =  (1.0/ 2.0)*WMNEF[  "mnab"]* R(2)[  "efmn"];

                XEIAM[  "eiam"]  = -           WMNEF[  "inaf"]* R(2)[  "efmn"];

                XIJMN[  "ijmn"]  =  (1.0/ 2.0)*WMNEF[  "ijef"]* R(2)[  "efmn"];

                XEIMN[  "eimn"] +=               FME[    "if"]* R(2)[  "efmn"];
                XEIMN[  "eimn"] +=             WMNEJ[  "oifn"]* R(2)[  "efmo"];
                XEIMN[  "eimn"] +=  (1.0/ 2.0)*WAMEF[  "eifg"]* R(2)[  "fgmn"];

                XEFAM[  "efam"] -=               FME[    "na"]* R(2)[  "efnm"];
                XEFAM[  "efam"] +=             WAMEF[  "enag"]* R(2)[  "fgmn"];
                XEFAM[  "efam"] +=  (1.0/ 2.0)*WMNEJ[  "noam"]* R(2)[  "efno"];

                Xi(2)[  "ijab"] +=               FME[    "ia"]*  DIA[    "jb"];
                Xi(2)[  "ijab"] +=             WAMEF[  "ejab"]*  DIA[    "ie"];
                Xi(2)[  "ijab"] -=             WMNEJ[  "ijam"]*  DIA[    "mb"];
                Xi(2)[  "ijab"] -=             WMNEF[  "imab"]*  DIJ[    "jm"];
                Xi(2)[  "ijab"] +=             WMNEF[  "ijae"]*  DAB[    "eb"];
                Xi(2)[  "ijab"] +=               XEA[    "eb"]* L(2)[  "ijae"];
                Xi(2)[  "ijab"] -=               XIM[    "jm"]* L(2)[  "imab"];
                Xi(2)[  "ijab"] +=  (1.0/ 2.0)*XEFAB[  "efab"]* L(2)[  "ijef"];
                Xi(2)[  "ijab"] -=             XEIAM[  "eiam"]* L(2)[  "mjeb"];
                Xi(2)[  "ijab"] +=  (1.0/ 2.0)*XIJMN[  "ijmn"]* L(2)[  "mnab"];
                Xi(2)[  "ijab"] +=  (1.0/ 2.0)*WMNEF[  "ijef"]*GABCD[  "efab"];
                Xi(2)[  "ijab"] -=             WMNEF[  "mjeb"]*GAIBJ[  "eiam"];
                Xi(2)[  "ijab"] +=  (1.0/ 2.0)*WMNEF[  "mnab"]*GIJKL[  "ijmn"];
                Xi(2)[  "ijab"] -=             WMNEJ[  "ijem"]*GAIBC[  "emab"];
                Xi(2)[  "ijab"] -=             WAMEF[  "emab"]*GIJAK[  "ijem"];
                Xi(2)[  "ijab"] -=             WAMEF[  "fiae"]*GAIBC[  "ejbf"];
                Xi(2)[  "ijab"] -=             WMNEJ[  "niam"]*GIJAK[  "mjbn"];
                Xi(2)[  "ijab"] +=               XEM[    "em"]* L(3)["ijmabe"];
                Xi(2)[  "ijab"] +=  (1.0/ 2.0)*XEFAM[  "efam"]* L(3)["ijmebf"];
                Xi(2)[  "ijab"] -=  (1.0/ 2.0)*XEIMN[  "eimn"]* L(3)["mnjeab"];

                  DIA[    "ia"] +=              L(2)[  "imae"]* R(1)[    "em"];

                  XIA[    "ia"]  =             WMNEF[  "imae"]* R(1)[    "em"];

                  XEA[    "ea"] -=               FME[    "ma"]* R(1)[    "em"];
                  XEA[    "ea"] +=             WAMEF[  "emaf"]* R(1)[    "fm"];

                  XIM[    "im"] +=               FME[    "ie"]* R(1)[    "em"];
                  XIM[    "im"] +=             WMNEJ[  "niem"]* R(1)[    "en"];

                GAIBC[  "aibc"] -=              L(2)[  "mibc"]* R(1)[    "am"];

                GIJAK[  "ijak"] +=              L(2)[  "ijae"]* R(1)[    "ek"];

                GIJAB[  "ijab"]  =              L(3)["ijmabe"]* R(1)[    "em"];

                XEFAB[  "efab"] -=             WAMEF[  "emab"]* R(1)[    "fm"];

                XEIAM[  "eiam"] +=             WAMEF[  "eiaf"]* R(1)[    "fm"];
                XEIAM[  "eiam"] -=             WMNEJ[  "niam"]* R(1)[    "en"];

                XIJMN[  "ijmn"] +=             WMNEJ[  "ijen"]* R(1)[    "em"];

                Xi(3)["ijkabc"] +=             WMNEF[  "ijab"]*  DIA[    "kc"];
                Xi(3)["ijkabc"] +=               FME[    "ia"]*GIJAB[  "jkbc"];
                Xi(3)["ijkabc"] +=               XIA[    "ia"]* L(2)[  "jkbc"];
                Xi(3)["ijkabc"] -=             WMNEF[  "imab"]*GIJAK[  "kjcm"];
                Xi(3)["ijkabc"] +=             WMNEF[  "ijeb"]*GAIBC[  "ekac"];
                Xi(3)["ijkabc"] -=             WMNEJ[  "kjcm"]*GIJAB[  "imab"];
                Xi(3)["ijkabc"] +=             WAMEF[  "ekac"]*GIJAB[  "ijeb"];
                Xi(3)["ijkabc"] +=               XEA[    "ea"]* L(3)["ijkebc"];
                Xi(3)["ijkabc"] -=               XIM[    "im"]* L(3)["mjkabc"];
                Xi(3)["ijkabc"] +=  (1.0/ 2.0)*XEFAB[  "efab"]* L(3)["ijkefc"];
                Xi(3)["ijkabc"] -=             XEIAM[  "eiam"]* L(3)["mjkebc"];
                Xi(3)["ijkabc"] +=  (1.0/ 2.0)*XIJMN[  "ijmn"]* L(3)["mnkabc"];

                printf("X3: %.15f\n", (1.0/36.0)*scalar(Xi(3)*Xi(3)));
                printf("X1: %.15f\n", Xi(1)({0,1},{1,0}).norm(2));
                printf("X2: %.15f\n", Xi(2)({0,1},{1,0}).norm(2));
            }

            return true;
        }
};

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::XiEOMEECCSDT);
REGISTER_TASK(aquarius::cc::XiEOMEECCSDT<double>,"xieomeeccsdt","");
