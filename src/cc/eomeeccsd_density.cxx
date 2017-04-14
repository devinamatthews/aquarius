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
class EOMEECCSDDensity : public Task
{
    public:
        EOMEECCSDDensity(const string& name, Config& config)
        : Task(name, config)
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("moints", "H"));
            reqs.push_back(Requirement("ccsd.T", "T"));
            reqs.push_back(Requirement("eomeeccsd.L", "L"));
            reqs.push_back(Requirement("eomeeccsd.R", "R"));
            reqs.push_back(Requirement("eomeeccsd.Zeta", "Z"));
            this->addProduct(Product("eomeeccsd.D", "D", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            auto& H = this->template get<TwoElectronOperator <U  >>("H");
            auto& T = this->template get<ExcitationOperator  <U,2>>("T");

            auto& Rs = this->template get<unique_vector<  ExcitationOperator<U,2>>>("R");
            auto& Ls = this->template get<unique_vector<DeexcitationOperator<U,2>>>("L");
            auto& Zs = this->template get<unique_vector<DeexcitationOperator<U,2>>>("Z");

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
                Z(1)[  "ia"] += R(0)*L(1)[  "ia"];
                Z(2)["ijab"] += R(0)*L(2)["ijab"];
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
            auto& XAMEF  = this->puttmp( "XAMEF", new SpinorbitalTensor<U>("X(am,ef)", arena, group, rep, {vrt,occ}, {1,1}, {2,0}));
            auto& XMNEJ  = this->puttmp( "XMNEJ", new SpinorbitalTensor<U>("X(mn,ej)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));
            auto& XAMEI  = this->puttmp( "XAMEI", new SpinorbitalTensor<U>("X(am,ei)", arena, group, rep, {vrt,occ}, {1,1}, {1,1}));
            auto& XMNIJ  = this->puttmp( "XMNIJ", new SpinorbitalTensor<U>("X(mn,ij)", arena, group, rep, {vrt,occ}, {0,2}, {0,2}));
            auto& XAMIJ  = this->puttmp( "XAMIJ", new SpinorbitalTensor<U>("X(am,ij)", arena, group, rep, {vrt,occ}, {1,1}, {0,2}));
            auto& GIJAK2 = this->puttmp("GIJAK2", new SpinorbitalTensor<U>("G(ij,ak)", arena, group, rep, {vrt,occ}, {0,2}, {1,1}));

            SpinorbitalTensor<U> Tau(T(2));
            Tau["abij"] += 0.5*T(1)["ai"]*T(1)["bj"];

              XAE[  "ae"]  = -0.5* L(2)["mnef"]*T(2)["afmn"];

              XMI[  "mi"]  =  0.5* L(2)["mnef"]*T(2)["efin"];

            XAMEF["amef"]  = -     L(2)["nmef"]*T(1)[  "an"];

            XMNEJ["mnej"]  =       L(2)["mnef"]*T(1)[  "fj"];

            XMNIJ["mnij"]  =  0.5* L(2)["mnef"]* Tau["efij"];

            XAMEI["amei"]  = -     L(2)["mnef"]*T(2)["afin"];
            XAMEI["amei"] -=      XMNEJ["nmei"]*T(1)[  "an"];

              DIA[  "ia"]  =       Z(1)[  "ia"];
              DIA[  "ia"] +=       L(2)["imae"]*R(1)[  "em"];

              DAB[  "ab"]  = -     L(1)[  "mb"]*R(1)[  "am"];
              DAB[  "ab"] -=  0.5* Z(2)["mnbe"]*T(2)["aemn"];
              DAB[  "ab"] -=  0.5* L(2)["mnbe"]*R(2)["aemn"];

              DIJ[  "ij"]  =        DIA[  "ie"]*T(1)[  "ej"];
              DIJ[  "ij"] +=       L(1)[  "ie"]*R(1)[  "ej"];
              DIJ[  "ij"] +=  0.5* Z(2)["imef"]*T(2)["efjm"];
              DIJ[  "ij"] +=  0.5* L(2)["imef"]*R(2)["efjm"];

              DAI[  "ai"]  =        DIA[  "me"]*T(2)["aeim"];
              DAI[  "ai"] +=       L(1)[  "me"]*R(2)["aeim"];
              DAI[  "ai"] -=        DIJ[  "mi"]*T(1)[  "am"];
              DAI[  "ai"] +=        DAB[  "ae"]*T(1)[  "ei"];
              DAI[  "ai"] -=        XMI[  "mi"]*R(1)[  "am"];
              DAI[  "ai"] +=        XAE[  "ae"]*R(1)[  "ei"];

              DAB[  "ab"] -=        DIA[  "mb"]*T(1)[  "am"];

              XAE[  "ae"] -=       L(1)[  "me"]*T(1)[  "am"];

              XAI[  "ai"]  =       L(1)[  "me"]*T(2)["aeim"];
              XAI[  "ai"] -=        XMI[  "mi"]*T(1)[  "am"];
              XAI[  "ai"] +=        XAE[  "ae"]*T(1)[  "ei"];

              XMI[  "mi"] +=       L(1)[  "me"]*T(1)[  "ei"];

            GIJAB["ijab"]  =       Z(2)["ijab"];

            GIJAK["ijak"]  =       Z(2)["ijae"]*T(1)[  "ek"];
            GIJAK2["ijak"] =       L(2)["ijae"]*R(1)[  "ek"];

            GAIBC["aibc"]  =      -Z(2)["mibc"]*T(1)[  "am"];
            GAIBC["aibc"] -=       L(2)["mibc"]*R(1)[  "am"];

            GIJKL["ijkl"]  =  0.5* Z(2)["ijef"]* Tau["efkl"];
            GIJKL["ijkl"] +=  0.5* L(2)["ijef"]*R(2)["efkl"];
            GIJKL["ijkl"] +=      XMNEJ["ijel"]*R(1)[  "ek"];

            GABCD["abcd"]  =  0.5* Z(2)["mncd"]* Tau["abmn"];
            GABCD["abcd"] +=  0.5* L(2)["mncd"]*R(2)["abmn"];
            GABCD["abcd"] -=      XAMEF["amcd"]*R(1)[  "bm"];

            GAIBJ["aibj"]  = -0.5* Z(2)["imbe"]*T(2)["aejm"];
            GAIBJ["aibj"] -=      GIJAK["mibj"]*T(1)[  "am"];
            GAIBJ["aibj"] -=      XMNEJ["mibj"]*R(1)[  "am"];
            GAIBJ["aibj"] +=      XAMEF["aibe"]*R(1)[  "ej"];

            GABCI["abci"]  = -      DIA[  "mc"]*T(2)["abmi"];
            GABCI["abci"] -=       L(1)[  "mc"]*R(2)["abmi"];
            GABCI["abci"] +=      GAIBC["amce"]*T(2)["beim"];
            GABCI["abci"] +=      XAMEF["amce"]*R(2)["beim"];
            GABCI["abci"] +=      GABCD["abce"]*T(1)[  "ei"];
            GABCI["abci"] += 0.5*GIJAK2["mnci"]* Tau["abmn"];
            GABCI["abci"] +=        DAB[  "ac"]*T(1)[  "bi"];
            GABCI["abci"] +=        XAE[  "ac"]*R(1)[  "bi"];

            GIJAK["ijak"] +=     GIJAK2["ijak"];

            GAIJK["aijk"]  =        DIA[  "ie"]*T(2)["aejk"];
            GAIJK["aijk"] +=       L(1)[  "ie"]*R(2)["aejk"];
            GAIJK["aijk"] +=      GIJAK["miek"]*T(2)["aejm"];
            GAIJK["aijk"] +=      XMNEJ["miek"]*R(2)["aejm"];
            GAIJK["aijk"] -=      GIJKL["mijk"]*T(1)[  "am"];
            GAIJK["aijk"] -=      XMNIJ["mijk"]*R(1)[  "am"];
            GAIJK["aijk"] +=        DIJ[  "ik"]*T(1)[  "aj"];
            GAIJK["aijk"] +=        XMI[  "ik"]*R(1)[  "aj"];

            XAMIJ["amij"]  =  0.5*XAMEF["amef"]* Tau["efij"];

            GABIJ["abij"]  =                     Tau["abij"];
            GABIJ["abij"] +=        DAI[  "ai"]*T(1)[  "bj"];
            GABIJ["abij"] +=        XAI[  "ai"]*R(1)[  "bj"];
            GABIJ["abij"] +=  0.5*GIJKL["mnij"]* Tau["abmn"];
            GABIJ["abij"] +=  0.5*XMNIJ["mnij"]*R(2)["abmn"];
            GABIJ["abij"] -=        DIJ[  "mi"]*T(2)["abmj"];
            GABIJ["abij"] +=        DAB[  "ae"]*T(2)["ebij"];
            GABIJ["abij"] -=        XMI[  "mi"]*R(2)["abmj"];
            GABIJ["abij"] +=        XAE[  "ae"]*R(2)["ebij"];
            GABIJ["abij"] -=      GAIBJ["amei"]*T(2)["ebmj"];
            GABIJ["abij"] -=      XAMEI["amei"]*R(2)["ebmj"];
            GABIJ["abij"] -=      XAMIJ["amij"]*R(1)[  "bm"];

            GAIBJ["aibj"] -=  0.5* Z(2)["imbe"]*T(2)["aejm"];
            GAIBJ["aibj"] -=       L(2)["imbe"]*R(2)["aejm"];
            GAIBJ["aibj"] -=        DIA[  "ib"]*T(1)[  "aj"];
            GAIBJ["aibj"] -=       L(1)[  "ib"]*R(1)[  "aj"];

              DAI[  "ai"] +=                    T(1)[  "ai"];

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

            U EIJAB = 0.25*scalar(GIJAB*H.getIJAB());
            U EIJAK =  0.5*scalar(GIJAK*H.getIJAK());
            U EAIBC =  0.5*scalar(GAIBC*H.getAIBC());
            U EABCD = 0.25*scalar(GABCD*H.getABCD());
            U EIJKL = 0.25*scalar(GIJKL*H.getIJKL());
            U EAIBJ =      scalar(GAIBJ*H.getAIBJ());
            U EAIJK =  0.5*scalar(GAIJK*H.getAIJK());
            U EABCI =  0.5*scalar(GABCI*H.getABCI());
            U EABIJ = 0.25*scalar(GABIJ*H.getABIJ());

            printf("|AI|: %20.15f\n", DAI({0,0},{0,0}).norm(2));
            printf("|IA|: %20.15f\n", DIA({0,0},{0,0}).norm(2));
            printf("|AB|: %20.15f\n", DAB({0,0},{0,0}).norm(2));
            printf("|IJ|: %20.15f\n", DIJ({0,0},{0,0}).norm(2));
            printf("|IJAK|: %20.15f\n", GIJAK({0,1},{0,1}).norm(2));
            printf("|AIJK|: %20.15f\n", GAIJK({0,1},{0,1}).norm(2));
            printf("|AIBC|: %20.15f\n", GAIBC({1,0},{1,0}).norm(2));
            printf("|ABCI|: %20.15f\n", GABCI({1,0},{1,0}).norm(2));

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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSDDensity);
REGISTER_TASK(aquarius::cc::EOMEECCSDDensity<double>,"eomeeccsd_density","");
