#include "util/global.hpp"

#include "time/time.hpp"
#include "task/task.hpp"
#include "operator/2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/deexcitationoperator.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;

namespace aquarius
{
namespace cc
{

template <typename U>
class EOMEECCSDTDensity : public Task
{
    public:
        EOMEECCSDTDensity(const string& name, Config& config)
        : Task(name, config)
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("moints", "H"));
            reqs.push_back(Requirement("ccsdt.T", "T"));
            reqs.push_back(Requirement("eomeeccsdt.L", "L"));
            reqs.push_back(Requirement("eomeeccsdt.R", "R"));
            reqs.push_back(Requirement("eomeeccsdt.Zeta", "Z"));
            this->addProduct(Product("eomeeccsdt.D", "D", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            auto& H = this->template get<TwoElectronOperator<U  >>("H");
            auto& T = this->template get<ExcitationOperator <U,3>>("T");

            auto& Rs = this->template get<unique_vector<  ExcitationOperator<U,3>>>("R");
            auto& Ls = this->template get<unique_vector<DeexcitationOperator<U,3>>>("L");
            auto& Zs = this->template get<unique_vector<DeexcitationOperator<U,3>>>("Z");

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            auto& D = this->put("D", new TwoElectronOperator<U>("D", arena, occ, vrt));

            auto& R = Rs.back();
            auto& L = Ls.back();
            auto& Z = Zs.back();

            auto& rep = R.getRepresentation();
            auto& group = rep.getPointGroup();

            if (rep)
            {
                Z(1)[    "ia"] += R(0)*L(1)[    "ia"];
                Z(2)[  "ijab"] += R(0)*L(2)[  "ijab"];
                Z(3)["ijkabc"] += R(0)*L(3)["ijkabc"];
            }

            auto& DIA = D.getIA();
            auto& DAI = D.getAI();
            auto& DAB = D.getAB();
            auto& DIJ = D.getIJ();

            auto& GIJAB = D.getIJAB();
            auto& GIJAK = D.getIJAK();
            auto& GAIBC = D.getAIBC();
            auto& GABCD = D.getABCD();
            auto& GIJKL = D.getIJKL();
            auto& GAIBJ = D.getAIBJ();
            auto& GAIJK = D.getAIJK();
            auto& GABCI = D.getABCI();
            auto& GABIJ = D.getABIJ();

            auto& XAE    = this->puttmp(   "XAE", new SpinorbitalTensor<U>(   "X(ae)", arena, group, rep, {vrt,occ}, {1,0}, {1,0}));
            auto& XMI    = this->puttmp(   "XMI", new SpinorbitalTensor<U>(   "X(mi)", arena, group, rep, {vrt,occ}, {0,1}, {0,1}));
            auto& XAI    = this->puttmp(   "XAI", new SpinorbitalTensor<U>(   "X(ai)", arena, group, rep, {vrt,occ}, {1,0}, {0,1}));
            auto& XABEF  = this->puttmp( "XABEF", new SpinorbitalTensor<U>("X(ab,ef)", arena, group, rep, {vrt,occ}, {2,0}, {2,0}));
            auto& XAMEF  = this->puttmp( "XAMEF", new SpinorbitalTensor<U>("X(am,ef)", arena, group, rep, {vrt,occ}, {1,1}, {2,0}));
            auto& XMNEJ  = this->puttmp( "XMNEJ", new SpinorbitalTensor<U>("X(mn,ej)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));
            auto& XAMEI  = this->puttmp( "XAMEI", new SpinorbitalTensor<U>("X(am,ei)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
            auto& XMNIJ  = this->puttmp( "XMNIJ", new SpinorbitalTensor<U>("X(mn,ij)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));
            auto& XABEJ  = this->puttmp( "XABEJ", new SpinorbitalTensor<U>("X(ab,ej)", arena, group, rep, {vrt,occ}, {2,0}, {1,1}));
            auto& XAMIJ  = this->puttmp( "XAMIJ", new SpinorbitalTensor<U>("X(am,ij)", arena, group, rep, {vrt,occ}, {1,1}, {0,2}));
            auto& GAIBC2 = this->puttmp("GAIBC2", new SpinorbitalTensor<U>("G(ai,bc)", arena, group, rep, {vrt,occ}, {1,1}, {2,0}));
            auto& GAIBJ2 = this->puttmp("GAIBJ2", new SpinorbitalTensor<U>("G(ai,bj)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
            auto& GIJAK2 = this->puttmp("GIJAK2", new SpinorbitalTensor<U>("G(ij,ak)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));

            printf("T1: %.15f\n", T(1)({1,0},{0,1}).norm(2));
            printf("T2: %.15f\n", T(2)({1,0},{0,1}).norm(2));
            printf("T3: %.15f\n", (1.0/36.0)*scalar(T(3)*T(3)));

            printf("R1: %.15f\n", R(1)({1,0},{0,1}).norm(2));
            printf("R2: %.15f\n", R(2)({1,0},{0,1}).norm(2));
            printf("R3: %.15f\n", (1.0/36.0)*scalar(R(3)*R(3)));

            printf("Z1: %.15f\n", Z(1)({0,1},{1,0}).norm(2));
            printf("Z2: %.15f\n", Z(2)({0,1},{1,0}).norm(2));
            printf("Z3: %.15f\n", (1.0/36.0)*scalar(Z(3)*Z(3)));

            printf("L1: %.15f\n", L(1)({0,1},{1,0}).norm(2));
            printf("L2: %.15f\n", L(2)({0,1},{1,0}).norm(2));
            printf("L3: %.15f\n", (1.0/36.0)*scalar(L(3)*L(3)));

            SpinorbitalTensor<U> Tau(T(2));
            Tau["abij"] += (1.0/ 2.0)*T(1)["ai"]*T(1)["bj"];

            SpinorbitalTensor<U> Taum(T(2));
            Taum["abij"] -= (1.0/ 2.0)*T(1)["ai"]*T(1)["bj"];

              DIA[  "ia"]  =              Z(1)[    "ia"];
              DIA[  "ia"] +=              L(2)[  "imae"]*R(1)[    "em"];
              DIA[  "ia"] +=  (1.0/ 4.0)* L(3)["imnaef"]*R(2)[  "efmn"];

            GIJAB["ijab"]  =              Z(2)[  "ijab"];
            GIJAB["ijab"] +=              L(3)["ijmabe"]*R(1)[    "em"];

              DAB[  "ab"]  = -            L(1)[    "mb"]*R(1)[    "am"];
              DAB[  "ab"] -=  (1.0/ 2.0)*GIJAB[  "mnbe"]*T(2)[  "aemn"];
              DAB[  "ab"] -=  (1.0/ 2.0)* L(2)[  "mnbe"]*R(2)[  "aemn"];
              DAB[  "ab"] -=  (1.0/12.0)* Z(3)["mnobef"]*T(3)["aefmno"];
              DAB[  "ab"] -=  (1.0/12.0)* L(3)["mnobef"]*R(3)["aefmno"];

              DIJ[  "ij"]  =               DIA[    "ie"]*T(1)[    "ej"];
              DIJ[  "ij"] +=              L(1)[    "ie"]*R(1)[    "ej"];
              DIJ[  "ij"] +=  (1.0/ 2.0)*GIJAB[  "imef"]*T(2)[  "efjm"];
              DIJ[  "ij"] +=  (1.0/ 2.0)* L(2)[  "imef"]*R(2)[  "efjm"];
              DIJ[  "ij"] +=  (1.0/12.0)* Z(3)["imnefg"]*T(3)["efgjmn"];
              DIJ[  "ij"] +=  (1.0/12.0)* L(3)["imnefg"]*R(3)["efgjmn"];

              XAE[  "ae"]  = -(1.0/ 2.0)* L(2)[  "mnef"]*T(2)[  "afmn"];
              XAE[  "ae"] -=  (1.0/12.0)* L(3)["mnoefg"]*T(3)["afgmno"];

              XMI[  "mi"]  =  (1.0/ 2.0)* L(2)[  "mnef"]*T(2)[  "efin"];
              XMI[  "mi"] +=  (1.0/12.0)* L(3)["mnoefg"]*T(3)["efgino"];

            GIJAK["ijak"]  =  (1.0/ 2.0)* Z(3)["ijmaef"]*T(2)[  "efkm"];
            GIJAK["ijak"] +=  (1.0/ 2.0)* L(3)["ijmaef"]*R(2)[  "efkm"];

            GAIBC["aibc"]  = -(1.0/ 2.0)* L(3)["minbce"]*R(2)[  "aemn"];

              DAI[  "ai"]  =               DIA[    "me"]*T(2)[  "aeim"];
              DAI[  "ai"] +=              L(1)[    "me"]*R(2)[  "aeim"];
              DAI[  "ai"] -=               DIJ[    "mi"]*T(1)[    "am"];
              DAI[  "ai"] +=               DAB[    "ae"]*T(1)[    "ei"];
              DAI[  "ai"] -=               XMI[    "mi"]*R(1)[    "am"];
              DAI[  "ai"] +=               XAE[    "ae"]*R(1)[    "ei"];
              DAI[  "ai"] -=  (1.0/ 2.0)*GIJAK[  "mnei"]*T(2)[  "eamn"];
              DAI[  "ai"] +=  (1.0/ 2.0)*GAIBC[  "amef"]*T(2)[  "efim"];
              DAI[  "ai"] +=  (1.0/ 4.0)*GIJAB[  "mnef"]*T(3)["aefimn"];
              DAI[  "ai"] +=  (1.0/ 4.0)* L(2)[  "mnef"]*R(3)["aefimn"];

              DAB[  "ab"] -=               DIA[    "mb"]*T(1)[    "am"];

            GAIBC["aibc"] -=  (1.0/ 2.0)* Z(3)["minbce"]*T(2)[  "aemn"];

            GIJKL["ijkl"]  =  (1.0/ 6.0)* Z(3)["ijmefg"]*T(3)["efgklm"];
            GIJKL["ijkl"] +=  (1.0/ 6.0)* L(3)["ijmefg"]*R(3)["efgklm"];
            GIJKL["ijkl"] +=             GIJAK[  "ijel"]*T(1)[    "ek"];

            GABCD["abcd"]  =  (1.0/ 6.0)* Z(3)["mnocde"]*T(3)["abemno"];
            GABCD["abcd"] +=  (1.0/ 6.0)* L(3)["mnocde"]*R(3)["abemno"];
            GABCD["abcd"] -=             GAIBC[  "amcd"]*T(1)[    "bm"];

            XMNEJ["mnej"]  =  (1.0/ 2.0)* L(3)["mnoefg"]*T(2)[  "fgjo"];
            XAMEF["amef"]  = -(1.0/ 2.0)* L(3)["nmoefg"]*T(2)[  "agno"];

            GIJAK["ijak"] +=             GIJAB[  "ijae"]*T(1)[    "ek"];

           GAIBJ2["aibj"]  = -           XMNEJ[  "mibj"]*R(1)[    "am"];
           GAIBJ2["aibj"] +=             XAMEF[  "aibe"]*R(1)[    "ej"];
           GAIBJ2["aibj"] -=  (1.0/ 4.0)* Z(3)["imnbef"]*T(3)["aefjmn"];
           GAIBJ2["aibj"] -=  (1.0/ 4.0)* L(3)["imnbef"]*R(3)["aefjmn"];

            GAIBJ["aibj"]  =            GAIBJ2[  "aibj"];
            GAIBJ["aibj"] -=  (1.0/ 2.0)*GIJAB[  "imbe"]*T(2)[  "aejm"];
            GAIBJ["aibj"] -=             GIJAK[  "mibj"]*T(1)[    "am"];
            GAIBJ["aibj"] +=             GAIBC[  "aibe"]*T(1)[    "ej"];
            GAIBJ["aibj"] +=             XMNEJ[  "mibj"]*R(1)[    "am"];
            GAIBJ["aibj"] -=             XAMEF[  "aibe"]*R(1)[    "ej"];

            GIJAK["ijak"] +=              L(2)[  "ijae"]*R(1)[    "ek"];
            GAIBC["aibc"] -=              L(2)[  "mibc"]*R(1)[    "am"];
            GAIBC["aibc"] -=             GIJAB[  "mibc"]*T(1)[    "am"];

            XABEF["abef"]  =  (1.0/ 6.0)* L(3)["mnoefg"]*T(3)["abgmno"];
            XABEF["abef"] -=             XAMEF[  "amef"]*T(1)[    "bm"];

            XAMEI["amei"]  = -(1.0/ 4.0)* L(3)["mnoefg"]*T(3)["afgino"];

            XMNIJ["mnij"]  =  (1.0/ 6.0)* L(3)["mnoefg"]*T(3)["efgijo"];
            XMNIJ["mnij"] +=             XMNEJ[  "mnej"]*T(1)[    "ei"];

              XAE[  "ae"] -=              L(1)[    "me"]*T(1)[    "am"];
              XMI[  "mi"] +=              L(1)[    "me"]*T(1)[    "ei"];

            XAMEF["amef"] -=              L(2)[  "nmef"]*T(1)[    "an"];
            XMNEJ["mnej"] +=              L(2)[  "mnef"]*T(1)[    "fj"];

            XABEJ["abci"]  =             XABEF[  "abce"]*T(1)[    "ei"];
            XABEJ["abci"] -=             XAMEI[  "amci"]*T(1)[    "bm"];
            XABEJ["abci"] +=             XAMEF[  "amce"]*T(2)[  "beim"];
            XABEJ["abci"] +=  (1.0/ 2.0)*XMNEJ[  "mnci"]* Tau[  "abmn"];
            XABEJ["abci"] -=  (1.0/ 2.0)* L(2)[  "mnce"]*T(3)["abemin"];

            XAMIJ["aijk"]  = -           XMNIJ[  "mijk"]*T(1)[    "am"];
            XAMIJ["aijk"] +=             XAMEI[  "aiek"]*T(1)[    "ej"];
            XAMIJ["aijk"] +=             XMNEJ[  "miek"]*T(2)[  "aejm"];
            XAMIJ["aijk"] +=  (1.0/ 2.0)*XAMEF[  "aief"]* Tau[  "efjk"];
            XAMIJ["aijk"] +=  (1.0/ 2.0)* L(2)[  "imef"]*T(3)["aefjkm"];

           GIJAK2["ijak"]  =              L(2)[  "ijae"]*T(1)[    "ek"];
            XAMEI["amei"] -=              L(2)[  "mnef"]*T(2)[  "afin"];
            XAMEI["amei"] -=            GIJAK2[  "nmei"]*T(1)[    "an"];

           GIJAK2["ijak"]  =             GIJAK[  "ijak"];
           GAIBC2["aibc"]  =             GAIBC[  "aibc"];

           GAIBC2["aibc"] +=              L(2)[  "mibc"]*R(1)[    "am"];
           GIJAK2["ijak"] -=              L(2)[  "ijae"]*R(1)[    "ek"];

            GABCI["abci"]  = -             DIA[    "mc"]*T(2)[  "abmi"];
            GABCI["abci"] -=              L(1)[    "mc"]*R(2)[  "abmi"];
            GABCI["abci"] +=            GAIBC2[  "amce"]*T(2)[  "beim"];
            GABCI["abci"] +=             XAMEF[  "amce"]*R(2)[  "beim"];
            GABCI["abci"] +=  (1.0/ 2.0)*GIJAK[  "mnci"]* Tau[  "abmn"];
            GABCI["abci"] +=  (1.0/ 2.0)*XMNEJ[  "mnci"]*R(2)[  "abmn"];
            GABCI["abci"] +=             GABCD[  "abce"]*T(1)[    "ei"];
            GABCI["abci"] +=             XABEF[  "abce"]*R(1)[    "ei"];
            GABCI["abci"] -=            GAIBJ2[  "amci"]*T(1)[    "bm"];
            GABCI["abci"] -=             XAMEI[  "amci"]*R(1)[    "bm"];
            GABCI["abci"] +=               DAB[    "ac"]*T(1)[    "bi"];
            GABCI["abci"] +=               XAE[    "ac"]*R(1)[    "bi"];
            GABCI["abci"] -=  (1.0/ 2.0)*GIJAB[  "mnce"]*T(3)["abemin"];
            GABCI["abci"] -=  (1.0/ 2.0)* L(2)[  "mnce"]*R(3)["abemin"];

            GAIJK["aijk"]  =               DIA[    "ie"]*T(2)[  "aejk"];
            GAIJK["aijk"] +=              L(1)[    "ie"]*R(2)[  "aejk"];
            GAIJK["aijk"] +=            GIJAK2[  "miek"]*T(2)[  "aejm"];
            GAIJK["aijk"] +=             XMNEJ[  "miek"]*R(2)[  "aejm"];
            GAIJK["aijk"] +=  (1.0/ 2.0)*GAIBC[  "aief"]* Tau[  "efjk"];
            GAIJK["aijk"] +=  (1.0/ 2.0)*XAMEF[  "aief"]*R(2)[  "efjk"];
            GAIJK["aijk"] -=             GIJKL[  "mijk"]*T(1)[    "am"];
            GAIJK["aijk"] -=             XMNIJ[  "mijk"]*R(1)[    "am"];
            GAIJK["aijk"] +=            GAIBJ2[  "aiek"]*T(1)[    "ej"];
            GAIJK["aijk"] +=             XAMEI[  "aiek"]*R(1)[    "ej"];
            GAIJK["aijk"] +=               DIJ[    "ik"]*T(1)[    "aj"];
            GAIJK["aijk"] +=               XMI[    "ik"]*R(1)[    "aj"];
            GAIJK["aijk"] +=  (1.0/ 2.0)*GIJAB[  "imef"]*T(3)["aefjkm"];
            GAIJK["aijk"] +=  (1.0/ 2.0)* L(2)[  "imef"]*R(3)["aefjkm"];

            GIJKL["ijkl"] +=  (1.0/ 2.0)*GIJAB[  "ijef"]* Tau[  "efkl"];
            GIJKL["ijkl"] +=  (1.0/ 2.0)* L(2)[  "ijef"]*R(2)[  "efkl"];

            XMNIJ["mnij"] +=  (1.0/ 2.0)* L(2)[  "mnef"]* Tau[  "efij"];

            XAMEF["amef"] +=              L(2)[  "nmef"]*T(1)[    "an"];
            XMNEJ["mnej"] -=              L(2)[  "mnef"]*T(1)[    "fj"];

              XMI[  "mi"] -=              L(1)[    "me"]*T(1)[    "ei"];

              XAI[  "ai"]  =              L(1)[    "me"]*T(2)[  "aeim"];
              XAI[  "ai"] -=               XMI[    "mi"]*T(1)[    "am"];
              XAI[  "ai"] +=               XAE[    "ae"]*T(1)[    "ei"];
              XAI[  "ai"] -=  (1.0/ 2.0)*XMNEJ[  "mnei"]*T(2)[  "eamn"];
              XAI[  "ai"] +=  (1.0/ 4.0)* L(2)[  "mnef"]*T(3)["aefimn"];

              XMI[  "mi"] +=              L(1)[    "me"]*T(1)[    "ei"];

            XAMEI["amei"] -=             XMNEJ[  "nmei"]*T(1)[    "an"];
            XAMEI["amei"] +=             XAMEF[  "amef"]*T(1)[    "fi"];

            GIJAK["mnej"] -=             GIJAB[  "mnef"]*T(1)[    "fj"];
            GIJAK["mnej"] -=              L(2)[  "mnef"]*R(1)[    "fj"];
            GAIBC["amef"] +=             GIJAB[  "nmef"]*T(1)[    "an"];
            GAIBC["amef"] +=              L(2)[  "nmef"]*R(1)[    "an"];

            GABIJ["abij"]  =               Tau[  "abij"];
            GABIJ["abij"] +=               DAI[    "ai"]*T(1)[    "bj"];
            GABIJ["abij"] +=               XAI[    "ai"]*R(1)[    "bj"];
            GABIJ["abij"] -=               DIJ[    "mi"]*T(2)[  "abmj"];
            GABIJ["abij"] +=               DAB[    "ae"]*T(2)[  "ebij"];
            GABIJ["abij"] +=  (1.0/ 2.0)*GABCD[  "abef"]* Tau[  "efij"];
            GABIJ["abij"] -=             GAIBJ[  "amei"]*T(2)[  "ebmj"];
            GABIJ["abij"] +=  (1.0/ 2.0)*GIJKL[  "mnij"]* Tau[  "abmn"];
            GABIJ["abij"] -=               XMI[    "mi"]*R(2)[  "abmj"];
            GABIJ["abij"] +=               XAE[    "ae"]*R(2)[  "ebij"];
            GABIJ["abij"] +=  (1.0/ 2.0)*XABEF[  "abef"]*R(2)[  "efij"];
            GABIJ["abij"] -=             XAMEI[  "amei"]*R(2)[  "ebmj"];
            GABIJ["abij"] +=  (1.0/ 2.0)*XMNIJ[  "mnij"]*R(2)[  "abmn"];
            GABIJ["abij"] +=             XABEJ[  "abej"]*R(1)[    "ei"];
            GABIJ["abij"] -=             XAMIJ[  "amij"]*R(1)[    "bm"];
            GABIJ["abij"] +=               DIA[    "me"]*T(3)["abeijm"];
            GABIJ["abij"] +=              L(1)[    "me"]*R(3)["abeijm"];
            GABIJ["abij"] +=  (1.0/ 2.0)*GAIBC[  "amef"]*T(3)["ebfijm"];
            GABIJ["abij"] -=  (1.0/ 2.0)*GIJAK[  "mnej"]*T(3)["abeinm"];
            GABIJ["abij"] +=  (1.0/ 2.0)*XAMEF[  "amef"]*R(3)["ebfijm"];
            GABIJ["abij"] -=  (1.0/ 2.0)*XMNEJ[  "mnej"]*R(3)["abeinm"];

            GIJAK["mnej"] +=             GIJAB[  "mnef"]*T(1)[    "fj"];
            GIJAK["mnej"] +=              L(2)[  "mnef"]*R(1)[    "fj"];
            GAIBC["amef"] -=             GIJAB[  "nmef"]*T(1)[    "an"];
            GAIBC["amef"] -=              L(2)[  "nmef"]*R(1)[    "an"];

            XAMEF["amef"] -=              L(2)[  "nmef"]*T(1)[    "an"];
            XMNEJ["mnej"] +=              L(2)[  "mnef"]*T(1)[    "fj"];

            GIJKL["ijkl"] +=             XMNEJ[  "ijel"]*R(1)[    "ek"];

            GABCD["abcd"] +=  (1.0/ 2.0)*GIJAB[  "mncd"]* Tau[  "abmn"];
            GABCD["abcd"] +=  (1.0/ 2.0)* L(2)[  "mncd"]*R(2)[  "abmn"];
            GABCD["abcd"] -=             XAMEF[  "amcd"]*R(1)[    "bm"];

            GAIBJ["aibj"] -=  (1.0/ 2.0)*GIJAB[  "imbe"]*T(2)[  "aejm"];
            GAIBJ["aibj"] -=             XMNEJ[  "mibj"]*R(1)[    "am"];
            GAIBJ["aibj"] +=             XAMEF[  "aibe"]*R(1)[    "ej"];
            GAIBJ["aibj"] -=              L(2)[  "imbe"]*R(2)[  "aejm"];
            GAIBJ["aibj"] -=               DIA[    "ib"]*T(1)[    "aj"];
            GAIBJ["aibj"] -=              L(1)[    "ib"]*R(1)[    "aj"];

              DAI[  "ai"] +=              T(1)[    "ai"];

              DAB *= -1;
              DIJ *= -1;
            GIJAK *= -1;
            GAIBC *= -1;
            GAIJK *= -1;
            GABCI *= -1;

            U EIA = scalar(DIA*H.getIA());
            U EAI = scalar(DAI*H.getAI());
            U EAB = scalar(DAB*H.getAB());
            U EIJ = scalar(DIJ*H.getIJ());

            U EIJAB = (1.0/4.0)*scalar(GIJAB*H.getIJAB());
            U EIJAK = (1.0/2.0)*scalar(GIJAK*H.getIJAK());
            U EAIBC = (1.0/2.0)*scalar(GAIBC*H.getAIBC());
            U EABCD = (1.0/4.0)*scalar(GABCD*H.getABCD());
            U EIJKL = (1.0/4.0)*scalar(GIJKL*H.getIJKL());
            U EAIBJ =           scalar(GAIBJ*H.getAIBJ());
            U EAIJK = (1.0/2.0)*scalar(GAIJK*H.getAIJK());
            U EABCI = (1.0/2.0)*scalar(GABCI*H.getABCI());
            U EABIJ = (1.0/4.0)*scalar(GABIJ*H.getABIJ());

            printf("|AI|: %20.15f\n", DAI({0,0},{0,0}).norm(2));
            printf("|IA|: %20.15f\n", DIA({0,0},{0,0}).norm(2));
            printf("|AB|: %20.15f\n", DAB({0,0},{0,0}).norm(2));
            printf("|IJ|: %20.15f\n", DIJ({0,0},{0,0}).norm(2));
            printf("|IJAK|: %20.15f\n", GIJAK({0,1},{0,1}).norm(2));
            printf("|AIJK|: %20.15f\n", GAIJK({0,1},{0,1}).norm(2));
            printf("|AIBC|: %20.15f\n", GAIBC({1,0},{1,0}).norm(2));
            printf("|ABCI|: %20.15f\n", GABCI({1,0},{1,0}).norm(2));
            printf("|AIBJ|: %20.15f\n", GAIBJ({1,0},{1,0}).norm(2));
            printf("|AIJB|: %20.15f\n", GAIBJ({1,0},{0,1}).norm(2));
            printf("|ABIJ|: %20.15f\n", GABIJ({1,0},{0,1}).norm(2));
            printf("|IJAB|: %20.15f\n", GIJAB({0,1},{1,0}).norm(2));
            printf("|ABCD|: %20.15f\n", GABCD({1,0},{1,0}).norm(2));
            printf("|IJKL|: %20.15f\n", GIJKL({0,1},{0,1}).norm(2));

            printf("\n");
            printf("EAB: %20.15f\n", EAB);
            printf("EIJ: %20.15f\n", EIJ);
            printf("EABCD: %20.15f\n", EABCD);
            printf("EIJKL: %20.15f\n", EIJKL);
            printf("EAIBJ: %20.15f\n", EAIBJ);
            printf("EAIJK: %20.15f\n", EAIJK+EIJAK);
            printf("EABCI: %20.15f\n", EABCI+EAIBC);
            printf("EABIJ: %20.15f\n", EABIJ+EIJAB);

            U E = EIA+EAI+EAB+EIJ+
                  EIJAB+EABIJ+EABCI+EAIBC+EAIJK+EIJAK+
                  EABCD+EIJKL+EAIBJ;

            printf("\n");
            printf("E: %20.15f\n", E);

            return true;
        }
};

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSDTDensity);
REGISTER_TASK(aquarius::cc::EOMEECCSDTDensity<double>,"eomeeccsdt_density","");
