#include "mp3.hpp"

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::time;
using namespace aquarius::task;

namespace aquarius
{
namespace cc
{

template <typename U>
MP3<U>::MP3(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("moints", "H"));
    addProduct(Product("double", "mp2", reqs));
    addProduct(Product("double", "mp3", reqs));
    addProduct(Product("double", "energy", reqs));
    addProduct(Product("double", "S2", reqs));
    addProduct(Product("double", "multiplicity", reqs));
}

template <typename U>
bool MP3<U>::run(TaskDAG& dag, const Arena& arena)
{
    const TwoElectronOperator<U>& H = get<TwoElectronOperator<U> >("H");

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    puttmp("T", new ExcitationOperator<U,2>("T", arena, occ, vrt));
    puttmp("D", new Denominator<U>(H));
    puttmp("Z", new ExcitationOperator<U,2>("Z", arena, occ, vrt));

    ExcitationOperator<U,2>& T = gettmp<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    Z(0) = (U)0.0;
    T(0) = (U)0.0;
    T(1) = (U)0.0;
    T(2) = H.getABIJ();

    T.weight(D);

    double energy = 0.25*real(scalar(H.getABIJ()*T(2)));
    double mp2energy = energy;


    Logger::log(arena) << "MP2 energy = " << setprecision(15) << energy << endl;
    put("mp2", new U(energy));

    TwoElectronOperator<U>& H_ = get<TwoElectronOperator<U> >("H");

    TwoElectronOperator<U> Hnew("W", H_, TwoElectronOperator<U>::AB|
                                 TwoElectronOperator<U>::IJ|
                                 TwoElectronOperator<U>::IJKL|
                                 TwoElectronOperator<U>::AIBJ);

    SpinorbitalTensor<U>& FAE = Hnew.getAB();
    SpinorbitalTensor<U>& FMI = Hnew.getIJ();
    SpinorbitalTensor<U>& WMNEF = Hnew.getIJAB();
    SpinorbitalTensor<U>& WABEF = Hnew.getABCD();
    SpinorbitalTensor<U>& WMNIJ = Hnew.getIJKL();
    SpinorbitalTensor<U>& WAMEI = Hnew.getAIBJ();

    Z(2)["abij"] = WMNEF["ijab"];
    Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
    Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*WABEF["abef"]*T(2)["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*T(2)["abmn"];
    Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];

    Z.weight(D);
    T += Z;

    energy = 0.25*real(scalar(H.getABIJ()*T(2)));

    Logger::log(arena) << "MP3 energy = " << setprecision(15) << energy << endl;
    Logger::log(arena) << "MP3 correlation energy = " << setprecision(15) << energy - mp2energy << endl;
    put("mp3", new U(energy));

    put("energy", new U(energy));

    /*
    if (isUsed("S2") || isUsed("multiplicity"))
    {
        double s2 = getProjectedS2(occ, vrt, T(1), T(2));
        double mult = sqrt(4*s2+1);

        put("S2", new Scalar(arena, s2));
        put("multiplicity", new Scalar(arena, mult));
    }
    */

    return true;
}

#if 0
template <typename U>
void MP3<U>::iterate(const Arena& arena)
{
    TwoElectronOperator<U>& H_ = get<TwoElectronOperator<U> >("H");

    ExcitationOperator<U,2>& T = get<ExcitationOperator<U,2> >("T");
    Denominator<U>& D = gettmp<Denominator<U> >("D");
    ExcitationOperator<U,2>& Z = gettmp<ExcitationOperator<U,2> >("Z");

    TwoElectronOperator<U> H("W", H_, TwoElectronOperator<U>::AB|
                                 TwoElectronOperator<U>::IJ|
                                 TwoElectronOperator<U>::IJKL|
                                 TwoElectronOperator<U>::AIBJ);

    SpinorbitalTensor<U>& FAE = H.getAB();
    SpinorbitalTensor<U>& FMI = H.getIJ();
    SpinorbitalTensor<U>& WMNEF = H.getIJAB();
    SpinorbitalTensor<U>& WABEF = H.getABCD();
    SpinorbitalTensor<U>& WMNIJ = H.getIJKL();
    SpinorbitalTensor<U>& WAMEI = H.getAIBJ();

//    sched.set_max_partitions(1);
    /**************************************************************************
     *
     * Intermediates
     */
    // FMI["mi"] += 0.5*WMNEF["mnef"]*T(2)["efin"];


    // WMNIJ["mnij"] += 0.5*WMNEF["mnef"]*T(2)["efij"];
    // FAE["ae"] -= 0.5*WMNEF["mnef"]*T(2)["afmn"];
    // WAMEI["amei"] -= 0.5*WMNEF["mnef"]*T(2)["afin"];
    /*
     *************************************************************************/

    /**************************************************************************
     *
     * T(1)->T(2) and T(2)->T(2)
     */
    Z(2)["abij"] = WMNEF["ijab"];
    Z(2)["abij"] += FAE["af"]*T(2)["fbij"];
    Z(2)["abij"] -= FMI["ni"]*T(2)["abnj"];
    Z(2)["abij"] += 0.5*WABEF["abef"]*T(2)["efij"];
    Z(2)["abij"] += 0.5*WMNIJ["mnij"]*T(2)["abmn"];
    Z(2)["abij"] -= WAMEI["amei"]*T(2)["ebmj"];
    /*
     *************************************************************************/

    Z.weight(D);
    T += Z;

    energy = 0.25*real(scalar(H.getABIJ()*T(2)));

    conv = Z.norm(00);

    diis.extrapolate(T, Z);
}
#endif

}
}

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::MP3);
REGISTER_TASK(aquarius::cc::MP3<double>,"mp3");
