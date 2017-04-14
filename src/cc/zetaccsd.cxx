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
class ZetaCCSD : public Iterative<U>
{
    protected:
        DIIS<DeexcitationOperator<U,2>> diis;

    public:
        ZetaCCSD(const string& name, Config& config)
        : Iterative<U>(name, config), diis(config.get("diis"))
        {
            vector<Requirement> reqs;
            reqs.push_back(Requirement("ccsd.Hbar", "Hbar"));
            reqs.push_back(Requirement("ccsd.T", "T"));
            reqs.push_back(Requirement(config.get<string>("source")+".Xi", "Xi"));
            this->addProduct(Product("double", "convergence", reqs));
            this->addProduct(Product(config.get<string>("source")+".Zeta", "Zeta", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            auto& Zetas = this->put("Zeta", new unique_vector<DeexcitationOperator<U,2>>());
            auto& D = this->puttmp("D", new Denominator<U>(H));
            auto& Z = this->puttmp("Z", new DeexcitationOperator<U,2>("Z", arena, occ, vrt));

            this->puttmp("GIM", new SpinorbitalTensor<U>("G(im)", H.getIJ()));
            this->puttmp("GEA", new SpinorbitalTensor<U>("G(ea)", H.getAB()));

            auto& Xis = this->template get<unique_vector<DeexcitationOperator<U,2>>>("Xi");

            for (int root = 0;root < Xis.size();root++)
            {
                auto& Xi = Xis[root];
                Zetas.emplace_back("Zeta", arena, occ, vrt);
                auto& Zeta = Zetas.back();

                Z   (0)         =             0;
                Zeta(0)         =             0;
                Zeta(1)[  "ia"] = Xi(1)[  "ia"];
                Zeta(2)["ijab"] = Xi(2)["ijab"];

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

            auto& Zetas = this->template get<unique_vector<DeexcitationOperator<U,2>>>("Zeta");
            auto& Xis   = this->template get<unique_vector<DeexcitationOperator<U,2>>>("Xi");

            auto& Zeta = Zetas[Zetas.size()-1];
            auto& Xi   =   Xis[Zetas.size()-1];

            auto& T = this->template get   <ExcitationOperator  <U,2>>("T");
            auto& D = this->template gettmp<Denominator         <U  >>("D");
            auto& Z = this->template gettmp<DeexcitationOperator<U,2>>("Z");

            auto& GIM = this->template gettmp<SpinorbitalTensor<U>>("GIM");
            auto& GEA = this->template gettmp<SpinorbitalTensor<U>>("GEA");

            /***************************************************************************
             *
             * Intermediates for Zeta-CCSD
             */
            GIM["mn"]  =  0.5*T(2)["efno"]*Zeta(2)["moef"];
            GEA["ef"]  = -0.5*T(2)["egmn"]*Zeta(2)["mnfg"];
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
            Z(1)[  "ia"] -=     WMNEJ["inam"]*    GIM[  "mn"];
            Z(1)[  "ia"] -=     WAMEF["fiea"]*    GEA[  "ef"];

            Z(2)["ijab"]  =     Xi(2)["ijab"];
            Z(2)["ijab"] +=       FME[  "ia"]*Zeta(1)[  "jb"];
            Z(2)["ijab"] +=     WAMEF["ejab"]*Zeta(1)[  "ie"];
            Z(2)["ijab"] -=     WMNEJ["ijam"]*Zeta(1)[  "mb"];
            Z(2)["ijab"] +=       FAE[  "ea"]*Zeta(2)["ijeb"];
            Z(2)["ijab"] -=       FMI[  "im"]*Zeta(2)["mjab"];
            Z(2)["ijab"] += 0.5*WABEF["efab"]*Zeta(2)["ijef"];
            Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]*Zeta(2)["mnab"];
            Z(2)["ijab"] +=     WAMEI["eiam"]*Zeta(2)["mjbe"];
            Z(2)["ijab"] -=     WMNEF["mjab"]*    GIM[  "im"];
            Z(2)["ijab"] +=     WMNEF["ijeb"]*    GEA[  "ea"];
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::ZetaCCSD);
REGISTER_TASK(aquarius::cc::ZetaCCSD<double>,"zetaccsd",spec);
