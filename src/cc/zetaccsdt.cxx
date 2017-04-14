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
class ZetaCCSDT : public Iterative<U>
{
    protected:
        DIIS<DeexcitationOperator<U,3>> diis;

    public:
        ZetaCCSDT(const string& name, Config& config)
        : Iterative<U>(name, config), diis(config.get("diis"))
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("ccsdt.Hbar", "Hbar"));
            reqs.push_back(Requirement("ccsdt.T", "T"));
            reqs.push_back(Requirement(config.get<string>("source")+".Xi", "Xi"));
            this->addProduct(Product("double", "convergence", reqs));
            this->addProduct(Product(config.get<string>("source")+".Zeta", "Zeta", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            auto& Zetas = this->put("Zeta", new unique_vector<DeexcitationOperator<U,3>>());
            auto& D = this->puttmp("D", new Denominator<U>(H));
            auto& Z = this->puttmp("Z", new DeexcitationOperator<U,3>("Z", arena, occ, vrt));

            this->puttmp(  "DIJ", new SpinorbitalTensor<U>(   "D(ij)",   H.getIJ()));
            this->puttmp(  "DAB", new SpinorbitalTensor<U>(   "D(ab)",   H.getAB()));
            this->puttmp(  "DAI", new SpinorbitalTensor<U>(   "D(ai)",   H.getAI()));
            this->puttmp("GABCD", new SpinorbitalTensor<U>("G(ab,cd)", H.getABCD()));
            this->puttmp("GAIBJ", new SpinorbitalTensor<U>("G(ai,bj)", H.getAIBJ()));
            this->puttmp("GIJKL", new SpinorbitalTensor<U>("G(ij,kl)", H.getIJKL()));
            this->puttmp("GAIBC", new SpinorbitalTensor<U>("G(ai,bc)", H.getAIBC()));
            this->puttmp("GIJAK", new SpinorbitalTensor<U>("G(ij,ak)", H.getIJAK()));

            auto& Xis = this->template get<unique_vector<DeexcitationOperator<U,3>>>("Xi");

            for (int root = 0;root < Xis.size();root++)
            {
                auto& Xi = Xis[root];
                Zetas.emplace_back("Zeta", arena, occ, vrt);
                auto& Zeta = Zetas.back();

                Z   (0)           =               0;
                Zeta(0)           =               0;
                Zeta(1)[    "ia"] = Xi(1)[    "ia"];
                Zeta(2)[  "ijab"] = Xi(2)[  "ijab"];
                Zeta(3)["ijkabc"] = Xi(3)["ijkabc"];

                Zeta.weight(D);

                Iterative<U>::run(dag, arena);

                this->put("convergence", new U(this->conv()));
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

            auto& Zetas = this->template get<unique_vector<DeexcitationOperator<U,3>>>("Zeta");
            auto& Xis   = this->template get<unique_vector<DeexcitationOperator<U,3>>>("Xi");

            auto& Zeta = Zetas[Zetas.size()-1];
            auto& Xi   =   Xis[Zetas.size()-1];

            auto& T = this->template get   <ExcitationOperator  <U,3>>("T");
            auto& D = this->template gettmp<Denominator         <U  >>("D");
            auto& Z = this->template gettmp<DeexcitationOperator<U,3>>("Z");

            auto&   DIJ = this->template gettmp<SpinorbitalTensor<U>>(  "DIJ");
            auto&   DAB = this->template gettmp<SpinorbitalTensor<U>>(  "DAB");
            auto&   DAI = this->template gettmp<SpinorbitalTensor<U>>(  "DAI");
            auto& GABCD = this->template gettmp<SpinorbitalTensor<U>>("GABCD");
            auto& GAIBJ = this->template gettmp<SpinorbitalTensor<U>>("GAIBJ");
            auto& GIJKL = this->template gettmp<SpinorbitalTensor<U>>("GIJKL");
            auto& GAIBC = this->template gettmp<SpinorbitalTensor<U>>("GAIBC");
            auto& GIJAK = this->template gettmp<SpinorbitalTensor<U>>("GIJAK");

            /***************************************************************************
             *
             * Intermediates for Zeta-CCSD
             */
            DIJ["ij"]  =  0.5*T(2)["efjm"]*Zeta(2)["imef"];
            DAB["ab"]  = -0.5*T(2)["aemn"]*Zeta(2)["mnbe"];
            /*
             **************************************************************************/

            /***************************************************************************
             *
             * Zeta-CCSD iteration
             */
            Z(1)[  "ia"]  =     Xi(1)[  "ia"];
            Z(1)[  "ia"] +=       FAE[  "ea"]*Zeta(1)[  "ie"];
            Z(1)[  "ia"] -=       FMI[  "im"]*Zeta(1)[  "ma"];
            Z(1)[  "ia"] -=     WAMEI["eiam"]*Zeta(1)[  "me"];
            Z(1)[  "ia"] += 0.5*WABEJ["efam"]*Zeta(2)["imef"];
            Z(1)[  "ia"] -= 0.5*WAMIJ["eimn"]*Zeta(2)["mnea"];
            Z(1)[  "ia"] -=     WMNEJ["inam"]*    DIJ[  "mn"];
            Z(1)[  "ia"] -=     WAMEF["fiea"]*    DAB[  "ef"];

            Z(2)["ijab"]  =     Xi(2)["ijab"];
            Z(2)["ijab"] +=       FME[  "ia"]*Zeta(1)[  "jb"];
            Z(2)["ijab"] +=     WAMEF["ejab"]*Zeta(1)[  "ie"];
            Z(2)["ijab"] -=     WMNEJ["ijam"]*Zeta(1)[  "mb"];
            Z(2)["ijab"] +=       FAE[  "ea"]*Zeta(2)["ijeb"];
            Z(2)["ijab"] -=       FMI[  "im"]*Zeta(2)["mjab"];
            Z(2)["ijab"] += 0.5*WABEF["efab"]*Zeta(2)["ijef"];
            Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*Zeta(2)["mnab"];
            Z(2)["ijab"] +=     WAMEI["eiam"]*Zeta(2)["mjbe"];
            Z(2)["ijab"] -=     WMNEF["mjab"]*    DIJ[  "im"];
            Z(2)["ijab"] +=     WMNEF["ijeb"]*    DAB[  "ea"];
            /*
             **************************************************************************/

            /***************************************************************************
             *
             * Intermediates for Zeta-CCSDT
             */
              DIJ[  "ij"]  =  (1.0/12.0)*T(3)["efgjmn"]*Zeta(3)["imnefg"];
              DAB[  "ab"]  = -(1.0/12.0)*T(3)["aefmno"]*Zeta(3)["mnobef"];

            GABCD["abcd"]  =   (1.0/6.0)*T(3)["abemno"]*Zeta(3)["mnocde"];
            GAIBJ["aibj"]  =       -0.25*T(3)["aefjmn"]*Zeta(3)["imnbef"];
            GIJKL["ijkl"]  =   (1.0/6.0)*T(3)["efgklm"]*Zeta(3)["ijmefg"];

            GIJAK["ijak"]  =         0.5*T(2)[  "efkm"]*Zeta(3)["ijmaef"];
            GAIBC["aibc"]  =        -0.5*T(2)[  "aemn"]*Zeta(3)["minbce"];

              DAI[  "ai"]  =        0.25*T(3)["aefimn"]*Zeta(2)[  "mnef"];
              DAI[  "ai"] -=         0.5*T(2)[  "eamn"]*  GIJAK[  "mnei"];
            /*
             **************************************************************************/

            /***************************************************************************
             *
             * Zeta-CCSDT iteration
             */
            Z(1)[    "ia"] +=       FME[  "ie"]*    DAB[    "ea"];
            Z(1)[    "ia"] -=       FME[  "ma"]*    DIJ[    "im"];
            Z(1)[    "ia"] -=     WMNEJ["inam"]*    DIJ[    "mn"];
            Z(1)[    "ia"] -=     WAMEF["fiea"]*    DAB[    "ef"];
            Z(1)[    "ia"] +=     WMNEF["miea"]*    DAI[    "em"];
            Z(1)[    "ia"] -= 0.5*WABEF["efga"]*  GAIBC[  "gief"];
            Z(1)[    "ia"] +=     WAMEI["eifm"]*  GAIBC[  "fmea"];
            Z(1)[    "ia"] -=     WAMEI["eman"]*  GIJAK[  "inem"];
            Z(1)[    "ia"] += 0.5*WMNIJ["imno"]*  GIJAK[  "noam"];
            Z(1)[    "ia"] -= 0.5*WAMEF["gief"]*  GABCD[  "efga"];
            Z(1)[    "ia"] +=     WAMEF["fmea"]*  GAIBJ[  "eifm"];
            Z(1)[    "ia"] -=     WMNEJ["inem"]*  GAIBJ[  "eman"];
            Z(1)[    "ia"] += 0.5*WMNEJ["noam"]*  GIJKL[  "imno"];

            Z(2)[  "ijab"] -=     WMNEF["mjab"]*    DIJ[    "im"];
            Z(2)[  "ijab"] +=     WMNEF["ijeb"]*    DAB[    "ea"];
            Z(2)[  "ijab"] += 0.5*WMNEF["ijef"]*  GABCD[  "efab"];
            Z(2)[  "ijab"] +=     WMNEF["imea"]*  GAIBJ[  "ejbm"];
            Z(2)[  "ijab"] += 0.5*WMNEF["mnab"]*  GIJKL[  "ijmn"];
            Z(2)[  "ijab"] -=     WAMEF["fiae"]*  GAIBC[  "ejbf"];
            Z(2)[  "ijab"] -=     WMNEJ["ijem"]*  GAIBC[  "emab"];
            Z(2)[  "ijab"] -=     WAMEF["emab"]*  GIJAK[  "ijem"];
            Z(2)[  "ijab"] -=     WMNEJ["niam"]*  GIJAK[  "mjbn"];
            Z(2)[  "ijab"] += 0.5*WABEJ["efbm"]*Zeta(3)["ijmaef"];
            Z(2)[  "ijab"] -= 0.5*WAMIJ["ejnm"]*Zeta(3)["imnabe"];

            Z(3)["ijkabc"]  =     Xi(3)["ijkabc"];
            Z(3)["ijkabc"] +=     WMNEF["ijab"]*Zeta(1)[    "kc"];
            Z(3)["ijkabc"] +=       FME[  "ia"]*Zeta(2)[  "jkbc"];
            Z(3)["ijkabc"] +=     WAMEF["ekbc"]*Zeta(2)[  "ijae"];
            Z(3)["ijkabc"] -=     WMNEJ["ijam"]*Zeta(2)[  "mkbc"];
            Z(3)["ijkabc"] +=     WMNEF["ijae"]*  GAIBC[  "ekbc"];
            Z(3)["ijkabc"] -=     WMNEF["mkbc"]*  GIJAK[  "ijam"];
            Z(3)["ijkabc"] +=       FAE[  "ea"]*Zeta(3)["ijkebc"];
            Z(3)["ijkabc"] -=       FMI[  "im"]*Zeta(3)["mjkabc"];
            Z(3)["ijkabc"] += 0.5*WABEF["efab"]*Zeta(3)["ijkefc"];
            Z(3)["ijkabc"] += 0.5*WMNIJ["ijmn"]*Zeta(3)["mnkabc"];
            Z(3)["ijkabc"] +=     WAMEI["eiam"]*Zeta(3)["mjkbec"];
            /*
             **************************************************************************/

            Z.weight(D);
            Zeta += Z;

            this->energy() = 0.0;
            this->conv() = Z.norm(00);

            diis.extrapolate(Zeta, Z);
        }
};

}
}

static const char* spec = R"!(

source
    string,
convergence?
    double 1e-9,
max_iterations?
    int 50,
conv_type?
    enum { MAXE, RMSE, MAE },
diis?
{
    damping?
        double 0.0,
    start?
        int 1,
    order?
        int 5,
    jacobi?
        bool false
}

)!";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::ZetaCCSDT);
REGISTER_TASK(aquarius::cc::ZetaCCSDT<double>,"zetaccsdt",spec);
