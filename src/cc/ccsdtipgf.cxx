#include "ccsdtipgf.hpp"

using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::input;
using namespace aquarius::op;
using namespace aquarius::convergence;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace cc
{

template <typename U>
CCSDTIPGF<U>::CCSDTIPGF(const string& name, Config& config)
: Iterative<CU>(name, config), krylov_config(config.get("krylov"))
{
    vector<Requirement> reqs;
    reqs.emplace_back("ccsdt.T", "T");
    reqs.emplace_back("ccsdt.L", "L");
    reqs.emplace_back("ccsdt.Hbar", "Hbar");
    this->addProduct("ccsdt.ipgf", "gf", reqs);

    orbital = config.get<int>("orbital");
    double from = config.get<double>("omega_min");
    double to = config.get<double>("omega_max");
    int n = config.get<double>("npoint");
    double eta = config.get<double>("eta");

    double delta = (to-from)/max(1,n-1);
    for (int i = 0;i < n;i++)
    {
        omegas.emplace_back(from+delta*i, eta);
    }
}

template <typename U>
bool CCSDTIPGF<U>::run(TaskDAG& dag, const Arena& arena)
{
    auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const PointGroup& group = H.getABIJ().getGroup();
    int nirrep = group.getNumIrreps();

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    int nI = occ.nalpha[0];
    int ni = occ.nbeta[0];
    int nA = vrt.nalpha[0];
    int na = vrt.nbeta[0];

    auto& T = this->template get<ExcitationOperator  <U,3>>("T");
    auto& L = this->template get<DeexcitationOperator<U,3>>("L");

    bool isalpha = false;
    bool isvrt = false;
    if (orbital > 0)
    {
        isalpha = true;
        orbital--;
        if (orbital >= nI)
        {
            isvrt = true;
            orbital -= nI;
        }
    }
    else
    {
        orbital = -orbital-1;
        if (orbital >= ni)
        {
            isvrt = true;
            orbital -= ni;
        }
    }

    auto& Rr = this->puttmp("Rr", new ExcitationOperator  <U,2,3>("Rr", arena, occ, vrt, isalpha ? -1 : 1));
    auto& Ri = this->puttmp("Ri", new ExcitationOperator  <U,2,3>("Ri", arena, occ, vrt, isalpha ? -1 : 1));
    auto& Zr = this->puttmp("Zr", new ExcitationOperator  <U,2,3>("Zr", arena, occ, vrt, isalpha ? -1 : 1));
    auto& Zi = this->puttmp("Zi", new ExcitationOperator  <U,2,3>("Zi", arena, occ, vrt, isalpha ? -1 : 1));
    auto& b  = this->puttmp("b",  new ExcitationOperator  <U,2,3>("b",  arena, occ, vrt, isalpha ? -1 : 1));
    auto& e  = this->puttmp("e",  new DeexcitationOperator<U,2,3>("e",  arena, occ, vrt, isalpha ? 1 : -1));

    auto& XE   = this->puttmp("XE",   new SpinorbitalTensor<U>("X(e)",    arena, group, {vrt,occ}, {0,0}, {1,0}, isalpha ? -1 : 1));
    auto& XMIJ = this->puttmp("XMIJ", new SpinorbitalTensor<U>("X(m,ij)", arena, group, {vrt,occ}, {0,1}, {0,2}, isalpha ? -1 : 1));
    auto& XAEI = this->puttmp("XAEI", new SpinorbitalTensor<U>("X(a,ei)", arena, group, {vrt,occ}, {1,0}, {1,1}, isalpha ? -1 : 1));
    auto& XMEI = this->puttmp("XMEI", new SpinorbitalTensor<U>("X(m,ei)", arena, group, {vrt,occ}, {0,1}, {1,1}, isalpha ? -1 : 1));
    auto& XAEF = this->puttmp("XAEF", new SpinorbitalTensor<U>("X(a,ef)", arena, group, {vrt,occ}, {1,0}, {2,0}, isalpha ? -1 : 1));

    SpinorbitalTensor<U> Dij    ("D(ij)",      arena, group, {vrt,occ}, {0,1}, {0,1});
    SpinorbitalTensor<U> Gijak  ("G(ij,ak)",   arena, group, {vrt,occ}, {0,2}, {1,1});
    SpinorbitalTensor<U> Gijkabl("G(ijk,abl)", arena, group, {vrt,occ}, {0,3}, {2,1});

        Dij[    "ij"]  =            L(1)[    "ie"]*T(1)[    "ej"];
        Dij[    "ij"] += (1.0/ 2.0)*L(2)[  "imef"]*T(2)[  "efjm"];
        Dij[    "ij"] += (1.0/12.0)*L(3)["imnefg"]*T(3)["efgjmn"];

      Gijak[  "ijak"]  =            L(2)[  "ijae"]*T(1)[    "ek"];
      Gijak[  "ijak"] += (1.0/ 2.0)*L(3)["ijmaef"]*T(2)[  "efkm"];

    Gijkabl["ijkabl"]  =            L(3)["ijkabe"]*T(1)[    "el"];

    SpinorbitalTensor<U> ap ("ap"  , arena, group, {vrt,occ}, {0,0}, {isvrt, !isvrt}, isalpha ? -1 : 1);
    SpinorbitalTensor<U> apt("ap^t", arena, group, {vrt,occ}, {isvrt, !isvrt}, {0,0}, isalpha ? 1 : -1);

    vector<tkv_pair<U>> pairs{{orbital, 1}};

    CTFTensor<U>& tensor1 = ap({0,0}, {isvrt && isalpha, !isvrt && isalpha})({0});
    if (arena.rank == 0)
        tensor1.writeRemoteData(pairs);
    else
        tensor1.writeRemoteData();

    CTFTensor<U>& tensor2 = apt({isvrt && isalpha, !isvrt && isalpha}, {0,0})({0});
    if (arena.rank == 0)
        tensor2.writeRemoteData(pairs);
    else
        tensor2.writeRemoteData();

    if (isvrt)
    {
        /*
         *  ab...    abe...
         * b  (e) = t
         *  ijk...   ijk...
         */
        b(1)[    "i"] = T(1)[    "ei"]*ap["e"];
        b(2)[  "aij"] = T(2)[  "aeij"]*ap["e"];
        b(3)["abijk"] = T(3)["abeijk"]*ap["e"];

        /*
         *  ijk...   ijk...
         * e  (e) = l
         *  ab...    abe...
         */
        e(1)[    "i"] = L(1)[    "ie"]*apt["e"];
        e(2)[  "ija"] = L(2)[  "ijae"]*apt["e"];
        e(3)["ijkab"] = L(3)["ijkabe"]*apt["e"];
    }
    else
    {
        /*
         * b (m) = d
         *  i       im
         */
        b(1)["i"] = ap["i"];

        /*
         *  ijk...           ij...     ijk...
         * e  (m) = d  (1 + l     ) + G
         *  ab...    km      ab...     abm...
         */
        e(1)[    "i"]  =                   apt["i"];
        e(1)[    "i"] -=     Dij[    "im"]*apt["m"];
        e(2)[  "ija"]  =    L(1)[    "ia"]*apt["j"];
        e(2)[  "ija"] -=   Gijak[  "ijam"]*apt["m"];
        e(3)["ijkab"]  =    L(2)[  "ijab"]*apt["k"];
        e(3)["ijkab"] -= Gijkabl["ijkabm"]*apt["m"];
    }

    //printf("<E1|E1>: %.15f\n",            scalar(e(1)*e(1)));
    //printf("<E2|E2>: %.15f\n", (1.0/ 2.0)*scalar(e(2)*e(2)));
    //printf("<E3|E3>: %.15f\n", (1.0/12.0)*scalar(e(3)*e(3)));
    //printf("<E|E>: %.15f\n", scalar(e*e));
    //printf("<B|B>: %.15f\n", scalar(b*b));
    //printf("<E|B>: %.15f\n",            scalar(e(1)[    "m"]*b(1)[    "m"]) +
    //                         (1.0/ 2.0)*scalar(e(2)[  "mne"]*b(2)[  "emn"]) +
    //                         (1.0/12.0)*scalar(e(3)["mnoef"]*b(3)["efmno"]));

    auto& D = this->puttmp("D", new ComplexDenominator<U>(H));

    for (auto& o : omegas)
    {
        this->puttmp("krylov", new ComplexLinearKrylov<ExcitationOperator<U,2,3>>(krylov_config, b));
        omega.real(-o.real());
        omega.imag( o.imag());

        this->log(arena) << "Computing Green's function at " << fixed << setprecision(6) << o << endl;

        Rr = b;
        Ri = 0;
        D.weight(Rr, Ri, omega);
        U norm = sqrt(aquarius::abs(scalar(Rr*Rr)) +
                      aquarius::abs(scalar(Ri*Ri)));
        Rr /= norm;
        Ri /= norm;

        Iterative<CU>::run(dag, arena);
    }

    return true;
}

template <typename U>
void CCSDTIPGF<U>::iterate(const Arena& arena)
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

    auto& XE   = this->template gettmp<SpinorbitalTensor<U>>("XE");
    auto& XMIJ = this->template gettmp<SpinorbitalTensor<U>>("XMIJ");
    auto& XAEI = this->template gettmp<SpinorbitalTensor<U>>("XAEI");
    auto& XMEI = this->template gettmp<SpinorbitalTensor<U>>("XMEI");
    auto& XAEF = this->template gettmp<SpinorbitalTensor<U>>("XAEF");

    auto& D = this->template gettmp<ComplexDenominator<U>>("D");
    auto& krylov = this->template gettmp<ComplexLinearKrylov<ExcitationOperator<U,2,3>>>("krylov");

    auto& Rr = this->template gettmp<  ExcitationOperator<U,2,3>>("Rr");
    auto& Ri = this->template gettmp<  ExcitationOperator<U,2,3>>("Ri");
    auto& Zr = this->template gettmp<  ExcitationOperator<U,2,3>>("Zr");
    auto& Zi = this->template gettmp<  ExcitationOperator<U,2,3>>("Zi");
    auto& b  = this->template gettmp<  ExcitationOperator<U,2,3>>("b");
    auto& e  = this->template gettmp<DeexcitationOperator<U,2,3>>("e");

    //printf("<Rr|Rr>: %.15f\n", scalar(Rr*Rr));
    //printf("<Ri|Ri>: %.15f\n", scalar(Ri*Ri));
    //printf("<Rr|Ri>: %.15f\n", scalar(Rr*Ri));

    //printf("<B|Rr>: %.15f\n", scalar(b*Rr));
    //printf("<B|Ri>: %.15f\n", scalar(b*Ri));

    for (int ri: {0,1})
    {
        ExcitationOperator<U,2,3>& R = {ri == 0 ? Rr : Ri};
        ExcitationOperator<U,2,3>& Z = {ri == 0 ? Zr : Zi};

          XE[    "e"]  = -0.5*WMNEF["mnfe"]*R(2)[   "fmn"];

        XMIJ[  "mij"]  =     -WMNIJ["mnij"]*R(1)[     "n"];
        XMIJ[  "mij"] +=      WMNEJ["nmei"]*R(2)[   "enj"];
        XMIJ[  "mij"] +=  0.5*WMNEF["mnef"]*R(3)[ "efinj"];

        XAEI[  "aei"]  =     -WAMEI["amei"]*R(1)[     "m"];
        XAEI[  "aei"] +=      WAMEF["amef"]*R(2)[   "fmi"];
        XAEI[  "aei"] +=  0.5*WMNEJ["mnei"]*R(2)[   "amn"];
        XAEI[  "aei"] -=  0.5*WMNEF["mnef"]*R(3)[ "afmni"];

        XMEI[  "mei"]  =     -WMNEJ["mnei"]*R(1)[     "n"];
        XMEI[  "mei"] +=      WMNEF["mnef"]*R(2)[   "fni"];

        XAEF[  "aef"]  =     -WAMEF["amef"]*R(1)[     "m"];
        XAEF[  "aef"] +=  0.5*WMNEF["mnef"]*R(2)[   "amn"];

        Z(1)[    "i"]  =       -FMI[  "mi"]*R(1)[     "m"];
        Z(1)[    "i"] +=        FME[  "me"]*R(2)[   "emi"];
        Z(1)[    "i"] -=  0.5*WMNEJ["mnei"]*R(2)[   "emn"];
        Z(1)[    "i"] += 0.25*WMNEF["mnef"]*R(3)[ "efmni"];

        Z(2)[  "aij"]  =     -WAMIJ["amij"]*R(1)[     "m"];
        Z(2)[  "aij"] +=        FAE[  "ae"]*R(2)[   "eij"];
        Z(2)[  "aij"] -=        FMI[  "mi"]*R(2)[   "amj"];
        Z(2)[  "aij"] +=  0.5*WMNIJ["mnij"]*R(2)[   "amn"];
        Z(2)[  "aij"] -=      WAMEI["amei"]*R(2)[   "emj"];
        Z(2)[  "aij"] +=         XE[   "e"]*T(2)[  "aeij"];
        Z(2)[  "aij"] +=        FME[  "me"]*R(3)[ "eamij"];
        Z(2)[  "aij"] +=  0.5*WAMEF["amef"]*R(3)[ "efimj"];
        Z(2)[  "aij"] -=  0.5*WMNEJ["mnej"]*R(3)[ "aeimn"];

        Z(3)["abijk"]  =      WABEJ["abej"]*R(2)[   "eik"];
        Z(3)["abijk"] -=      WAMIJ["amij"]*R(2)[   "bmk"];
        Z(3)["abijk"] -=       XMIJ[ "mik"]*T(2)[  "abmj"];
        Z(3)["abijk"] -=       XAEI[ "aei"]*T(2)[  "bejk"];
        Z(3)["abijk"] +=        FAE[  "ae"]*R(3)[ "ebijk"];
        Z(3)["abijk"] -=        FMI[  "mi"]*R(3)[ "abmjk"];
        Z(3)["abijk"] -=      WAMEI["amei"]*R(3)[ "ebmjk"];
        Z(3)["abijk"] +=  0.5*WABEF["abef"]*R(3)[ "efijk"];
        Z(3)["abijk"] +=  0.5*WMNIJ["mnij"]*R(3)[ "abmnk"];
        Z(3)["abijk"] +=         XE[   "e"]*T(3)["abeijk"];
        Z(3)["abijk"] +=       XMEI[ "mek"]*T(3)["abeijm"];
        Z(3)["abijk"] +=   0.5*XAEF[ "bef"]*T(3)["aefijk"];
    }

    //printf("<Z1|Z1>: %.15f\n", scalar(Zr(1)*Zr(1)));
    //printf("<Z2|Z2>: %.15f\n", 0.5*scalar(Zr(2)*Zr(2)));
    //printf("<Zr|Zr>: %.15f\n", scalar(Zr*Zr));
    //printf("<Zi|Zi>: %.15f\n", scalar(Zi*Zi));

    /*
     * Convert H*r to (H-w)*r
     */
    Zr -= omega.real()*Rr;
    Zr += omega.imag()*Ri;
    Zi -= omega.real()*Ri;
    Zi -= omega.imag()*Rr;
    //Zr *= -1;
    //Zi *= -1;

    //printf("<Ur|Ur>: %.15f\n", scalar(Zr*Zr));
    //printf("<Ui|Ui>: %.15f\n", scalar(Zi*Zi));

    krylov.extrapolate(Rr, Ri, Zr, Zi, D, omega);

    this->conv() = max(Zr.norm(00), Zi.norm(00));

    krylov.getSolution(Zr, Zi);

    this->energy() = CU(           scalar(e(1)[    "m"]*Zr(1)[    "m"]) +
                        (1.0/ 2.0)*scalar(e(2)[  "mne"]*Zr(2)[  "emn"]) +
                        (1.0/12.0)*scalar(e(3)["mnoef"]*Zr(3)["efmno"]),
                                   scalar(e(1)[    "m"]*Zi(1)[    "m"]) +
                        (1.0/ 2.0)*scalar(e(2)[  "mne"]*Zi(2)[  "emn"]) +
                        (1.0/12.0)*scalar(e(3)["mnoef"]*Zi(3)["efmno"]));
}

}
}

static const char* spec = R"(

orbital int,
npoint int,
omega_min double,
omega_max double,
eta double,
convergence?
    double 1e-9,
max_iterations?
    int 150,
conv_type?
    enum { MAXE, RMSE, MAE },
krylov?
{
    order?
            int 10,
    compaction?
            enum { discrete, continuous },
}

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::CCSDTIPGF);
REGISTER_TASK(aquarius::cc::CCSDTIPGF<double>, "ccsdtipgf",spec);
