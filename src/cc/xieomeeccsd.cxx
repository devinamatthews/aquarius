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
class XiEOMEECCSD : public Task
{
    public:
        XiEOMEECCSD(const string& name, Config& config)
        : Task(name, config)
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
            reqs.push_back(Requirement("eomeeccsd.R", "R"));
            reqs.push_back(Requirement("eomeeccsd.L", "L"));
            this->addProduct(Product("eomeeccsd.Xi", "Xi", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;
            auto& group = occ.group;

            auto& Xis = this->put("Xi", new unique_vector<DeexcitationOperator<U,2>>());

            auto& Rs = this->template get<unique_vector<ExcitationOperator<U,2>>>("R");
            auto& Ls = this->template get<unique_vector<ExcitationOperator<U,2>>>("L");

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
                auto& GIJAK = this->puttmp("GIJAK", new SpinorbitalTensor<U>("G(ij,ak)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));

                auto&   XIA = this->puttmp(  "XIA", new SpinorbitalTensor<U>(   "X(ia)", arena, group, rep, {vrt,occ}, {0,1}, {1,0}));
                auto&   XEM = this->puttmp(  "XEM", new SpinorbitalTensor<U>(   "X(em)", arena, group, rep, {vrt,occ}, {1,0}, {0,1}));
                auto&   XEA = this->puttmp(  "XEA", new SpinorbitalTensor<U>(   "X(ea)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
                auto&   XIM = this->puttmp(  "XIM", new SpinorbitalTensor<U>(   "X(im)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
                auto& XEIAM = this->puttmp("XEIAM", new SpinorbitalTensor<U>("X(ei,am)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
                auto& XIJMN = this->puttmp("XIJMN", new SpinorbitalTensor<U>("X(ij,mn)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));
                auto& XEIMN = this->puttmp("XEIMN", new SpinorbitalTensor<U>("X(ei,mn)", arena, group, rep, {vrt,occ}, {1,1}, {0,2}));
                auto& XEFAM = this->puttmp("XEFAM", new SpinorbitalTensor<U>("X(ef,am)", arena, group, rep, {vrt,occ}, {2,0}, {1,1}));

                U LR2 = U(1) - scalar(L(1)["ia"]*R(1)["ai"]);

                Xi(0)         =                 0;
                Xi(1)[  "ia"] = LR2*  FME[  "ia"];
                Xi(2)["ijab"] =     WMNEF["ijab"];

                if (rep)
                {
                    U FR = scalar(FME["ia"]*R(1)["ai"]);
                    U WR = 0.25*scalar(WMNEF["ijab"]*R(2)["abij"]);

                    Xi(1)[  "ia"] +=     WR *L(1)[  "ia"];
                    Xi(2)["ijab"] += (FR+WR)*L(2)["ijab"];
                }

                  XIA[  "ia"]  =      WMNEF["imae"]* R(1)[  "em"];

                  XEM[  "em"]  =        FME[  "nf"]* R(2)["efmn"];
                  XEM[  "em"] +=  0.5*WAMEF["enfg"]* R(2)["fgmn"];
                  XEM[  "em"] -=  0.5*WMNEJ["nofm"]* R(2)["feno"];

                  XEA[  "ea"]  = -0.5*WMNEF["mnaf"]* R(2)["efmn"];

                  XIM[  "im"]  =  0.5*WMNEF["inef"]* R(2)["efmn"];

                  DIA[  "ia"]  =       L(2)["imae"]* R(1)[  "em"];

                  DAI[  "ai"]  =       L(1)[  "me"]* R(2)["aeim"];

                  DIJ[  "ij"]  =  0.5* L(2)["imef"]* R(2)["efjm"];

                  DAB[  "ab"]  = -0.5* L(2)["mnbe"]* R(2)["aemn"];

                GIJKL["ijkl"]  =  0.5* L(2)["ijef"]* R(2)["efkl"];

                GIJAK["ijak"]  =       L(2)["ijae"]* R(1)[  "ek"];

                XEIMN["eimn"]  =        FME[  "if"]* R(2)["efmn"];
                XEIMN["eimn"] +=      WMNEJ["oifn"]* R(2)["efmo"];
                XEIMN["eimn"] +=  0.5*WAMEF["eifg"]* R(2)["fgmn"];

                XEFAM["efam"]  = -      FME[  "na"]* R(2)["efnm"];
                XEFAM["efam"] +=      WAMEF["enag"]* R(2)["fgmn"];
                XEFAM["efam"] +=  0.5*WMNEJ["noam"]* R(2)["efno"];

                XEIAM["eiam"]  = -    WMNEF["inaf"]* R(2)["efmn"];
                XEIAM["eiam"] +=      WAMEF["eiaf"]* R(1)[  "fm"];
                XEIAM["eiam"] -=      WMNEJ["niam"]* R(1)[  "en"];

                XIJMN["ijmn"]  =  0.5*WMNEF["ijef"]* R(2)["efmn"];
                XIJMN["ijmn"] +=      WMNEJ["ijen"]* R(1)[  "em"];

                Xi(1)[  "ia"] +=        XEA[  "ea"]* L(1)[  "ie"];
                Xi(1)[  "ia"] -=        XIM[  "im"]* L(1)[  "ma"];
                Xi(1)[  "ia"] +=      WMNEF["imae"]*  DAI[  "em"];
                Xi(1)[  "ia"] -=      WMNEJ["iman"]*  DIJ[  "nm"];
                Xi(1)[  "ia"] -=      WAMEF["eifa"]*  DAB[  "fe"];
                Xi(1)[  "ia"] +=        XEM[  "em"]* L(2)["imae"];
                Xi(1)[  "ia"] +=  0.5*XEFAM["efam"]* L(2)["imef"];
                Xi(1)[  "ia"] -=  0.5*XEIMN["eimn"]* L(2)["mnea"];

                  XEA[  "ea"] -=        FME[  "ma"]* R(1)[  "em"];
                  XEA[  "ea"] +=      WAMEF["emaf"]* R(1)[  "fm"];

                  XIM[  "im"] +=        FME[  "ie"]* R(1)[  "em"];
                  XIM[  "im"] +=      WMNEJ["niem"]* R(1)[  "en"];

                  DIJ[  "ij"] +=       L(1)[  "ie"]* R(1)[  "ej"];

                  DAB[  "ab"] -=       L(1)[  "mb"]* R(1)[  "am"];

                Xi(2)["ijab"] +=        XIA[  "ia"]* L(1)[  "jb"];
                Xi(2)["ijab"] +=        FME[  "ia"]*  DIA[  "jb"];
                Xi(2)["ijab"] +=      WAMEF["ejab"]*  DIA[  "ie"];
                Xi(2)["ijab"] -=      WMNEJ["ijam"]*  DIA[  "mb"];
                Xi(2)["ijab"] -=      WMNEF["imab"]*  DIJ[  "jm"];
                Xi(2)["ijab"] +=      WMNEF["ijae"]*  DAB[  "eb"];
                Xi(2)["ijab"] +=        XEA[  "eb"]* L(2)["ijae"];
                Xi(2)["ijab"] -=        XIM[  "jm"]* L(2)["imab"];
                Xi(2)["ijab"] -=      XEIAM["eiam"]* L(2)["mjeb"];
                Xi(2)["ijab"] +=  0.5*XIJMN["ijmn"]* L(2)["mnab"];
                Xi(2)["ijab"] +=  0.5*WMNEF["mnab"]*GIJKL["ijmn"];
                Xi(2)["ijab"] -=      WAMEF["emab"]*GIJAK["ijem"];
            }

            return true;
        }
};

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::XiEOMEECCSD);
REGISTER_TASK(aquarius::cc::XiEOMEECCSD<double>,"xieomeeccsd","");
