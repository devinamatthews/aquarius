#include "ccsdipgf.hpp"

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
CCSDIPGF<U>::CCSDIPGF(const string& name, Config& config)
: Iterative<U>(name, config)
{
    vector<Requirement> reqs;
    reqs.emplace_back("ccsd.T", "T");
    reqs.emplace_back("ccsd.L", "L");
    reqs.emplace_back("ccsd.Hbar", "Hbar");
    this->addProduct("ccsd.ipgf", "gf", reqs);

    orbital = config.get<int>("orbital");
    double from = config.get<double>("omega.from");
    double to = config.get<double>("omega.to");
    int n = config.get<double>("omega.num_points");
    double eta = config.get<double>("eta");

    double delta = (to-from)/max(1,n-1);
    for (int i = 0;i < eta;i++)
    {
        omega.emplace_back(from+delta*i, eta);
    }
}

template <typename U>
bool CCSDIPGF<U>::run(TaskDAG& dag, const Arena& arena)
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

    this->puttmp("D", new Denominator<U>(H));

    auto& T = this->template get<ExcitationOperator  <U,2>>("T");
    auto& L = this->template get<DeexcitationOperator<U,2>>("L");

    if (orbital > 0)
    {
        isalpha = true;
        orbital--;
        if (orbital > nI)
        {
            isvrt = true;
            orbital -= nI;
        }
    }
    else
    {
        orbital = -orbital-1;
        if (orbital > ni)
        {
            isvrt = true;
            orbital -= ni;
        }
    }

    auto& Xr = this->puttmp("Xr", new ExcitationOperator  <U,1,2>("Xr", arena, occ, vrt, isalpha ? 1 : -1));
    auto& Xi = this->puttmp("Xi", new ExcitationOperator  <U,1,2>("Xi", arena, occ, vrt, isalpha ? 1 : -1));
    auto& b  = this->puttmp("b",  new ExcitationOperator  <U,1,2>("b",  arena, occ, vrt, isalpha ? 1 : -1));
    auto& e  = this->puttmp("e",  new DeexcitationOperator<U,1,2>("e",  arena, occ, vrt, isalpha ? 1 : -1));

    SpinorbitalTensor<U> Dij("D(ij)", arena, group, {vrt,occ}, {0,1}, {0,1});
    SpinorbitalTensor<U> Gijak("G(ij,ak)", arena, group, {vrt,occ}, {0,2}, {1,1});

    Dij["ij"]  =     L(1)["ie"  ]*T(1)["ej"  ];
    Dij["ij"] += 0.5*L(2)["imef"]*T(2)["efjm"];

    Gijak["ijak"] = L(2)["ijae"]*T(1)["ek"];

    SpinorbitalTensor<U> ap ("ap"  , arena, group, {vrt,occ}, {0,0}, {isvrt ? 1 : 0, isvrt ? 0 : 1}, isalpha ? 1 : -1);
    SpinorbitalTensor<U> apt("ap^t", arena, group, {vrt,occ}, {isvrt ? 1 : 0, isvrt ? 0 : 1}, {0,0}, isalpha ? 1 : -1);

    vector<tkv_pair<U>> pairs{{orbital, 1}};

    CTFTensor<U>& tensor1 = ap({0,0}, {isvrt ? (isalpha ? 1 : 0) : 0, isvrt ? 0 : (isalpha ? 1 : 0)})({0});
    if (arena.rank == 0)
        tensor1.writeRemoteData(pairs);
    else
        tensor1.writeRemoteData();

    CTFTensor<U>& tensor2 = apt({isvrt ? (isalpha ? 1 : 0) : 0, isvrt ? 0 : (isalpha ? 1 : 0)}, {0,0})({0});
    if (arena.rank == 0)
        tensor2.writeRemoteData(pairs);
    else
        tensor2.writeRemoteData();

    if (isvrt)
    {
        /*
         *  ab...    ab...e
         * b  (e) = t
         *  ijk...   ijk...
         */
        b(1)["i"  ] = T(1)["ei"  ]*ap["e"];
        b(2)["aij"] = T(2)["aeij"]*ap["e"];

        /*
         *  ijk...   ijk...
         * e  (e) = l
         *  ab...    ab...e
         */
        e(1)["i"  ] = L(1)["ie"  ]*apt["e"];
        e(2)["ija"] = L(2)["ijae"]*apt["e"];
    }
    else
    {
        /*
         * b (m) = d
         *  i       im
         */
        b(1)["i"] = ap["i"];

        /*
         *  ijk...   ijk...
         * e  (m) = G
         *  ab...    ab...m
         */
        e(1)["i"  ] = Dij  ["im"  ]*apt["m"];
        e(2)["ija"] = Gijak["ijam"]*apt["m"];
    }

    return true;
}

template <typename U>
void CCSDIPGF<U>::iterate(const Arena& arena)
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

    auto& T = this->template get<ExcitationOperator<U,2>>("T");

    auto& XMI = this->template gettmp<SpinorbitalTensor<U>>("XMI");
    auto& XAE = this->template gettmp<SpinorbitalTensor<U>>("XAE");

    auto& D = this->template gettmp<Denominator<U>>("D");
    auto& davidson = this->template gettmp<Davidson<ExcitationOperator<U,2>>>("Davidson");

    auto& Rs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("R");
    auto& Zs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("Z");
    auto& Vs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("V");

    for (int root = 0;root < this->nsolution();root++)
    {
        ExcitationOperator<U,2>& R = Rs[root];
        ExcitationOperator<U,2>& Z = Zs[root];
        Z = 0;

        for (auto& V : Vs)
        {
            R -= scalar(conj(R)*V)*V;
        }
        R /= sqrt(aquarius::abs(scalar(conj(R)*R)));

         XMI[  "mi"]  =     WMNEJ["nmei"]*R(1)[  "en"];
         XMI[  "mi"] += 0.5*WMNEF["mnef"]*R(2)["efin"];
         XAE[  "ae"]  =     WAMEF["amef"]*R(1)[  "fm"];
         XAE[  "ae"] -= 0.5*WMNEF["mnef"]*R(2)["afmn"];

        Z(1)[  "ai"] +=       FAE[  "ae"]*R(1)[  "ei"];
        Z(1)[  "ai"] -=       FMI[  "mi"]*R(1)[  "am"];
        Z(1)[  "ai"] -=     WAMEI["amei"]*R(1)[  "em"];
        Z(1)[  "ai"] +=       FME[  "me"]*R(2)["aeim"];
        Z(1)[  "ai"] += 0.5*WAMEF["amef"]*R(2)["efim"];
        Z(1)[  "ai"] -= 0.5*WMNEJ["mnei"]*R(2)["eamn"];

        Z(2)["abij"] +=     WABEJ["abej"]*R(1)[  "ei"];
        Z(2)["abij"] -=     WAMIJ["amij"]*R(1)[  "bm"];
        Z(2)["abij"] +=       FAE[  "ae"]*R(2)["ebij"];
        Z(2)["abij"] -=       FMI[  "mi"]*R(2)["abmj"];
        Z(2)["abij"] +=       XAE[  "ae"]*T(2)["ebij"];
        Z(2)["abij"] -=       XMI[  "mi"]*T(2)["abmj"];
        Z(2)["abij"] += 0.5*WMNIJ["mnij"]*R(2)["abmn"];
        Z(2)["abij"] += 0.5*WABEF["abef"]*R(2)["efij"];
        Z(2)["abij"] -=     WAMEI["amei"]*R(2)["ebmj"];

        //0.5*Z(1)({0,0},{0,0})[  "ai"] += 0.5*Z(1)({1,0},{0,1})[  "ai"];
        //    Z(1)({1,0},{0,1})[  "ai"]  =     Z(1)({0,0},{0,0})[  "ai"];
        //0.5*Z(2)({1,0},{0,1})["abij"] += 0.5*Z(2)({1,0},{0,1})["baji"];
        //0.5*Z(2)({0,0},{0,0})["abij"] += 0.5*Z(2)({2,0},{0,2})["abij"];
        //    Z(2)({2,0},{0,2})["abij"]  =     Z(2)({0,0},{0,0})["abij"];
    }

    vector<U> energies = davidson.extrapolate(Rs, Zs, D);

    for (int i = 0;i < this->nsolution();i++)
    {
        this->energy(i) = energies[i];
        this->conv(i) = Zs[i].norm(00);
    }
}

}
}

static const char* spec = R"(

multiroot?
    bool false,
nroot?
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
            int 500,
    jacobi?
            bool false
}

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSD);
REGISTER_TASK(aquarius::cc::EOMEECCSD<double>, "eomeeccsd",spec);
