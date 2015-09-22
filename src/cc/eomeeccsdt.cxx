#include "eomeeccsdt.hpp"

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
EOMEECCSDT<U>::EOMEECCSDT(const string& name, Config& config)
: Iterative<U>(name, config), davidson_config(config.get("davidson")),
  nroot(config.get<int>("nroot")), ntriplet(config.get<int>("ntriplet")), multiroot(config.get<bool>("multiroot"))
{
    vector<Requirement> reqs;
    reqs.emplace_back("ccsdt.T", "T");
    reqs.emplace_back("ccsdt.Hbar", "Hbar");
    reqs.emplace_back("tda.TDAevals", "TDAevals");
    reqs.emplace_back("tda.TDAevecs", "TDAevecs");
    this->addProduct("eomeeccsdt.energy", "energy", reqs);
    this->addProduct("eomeeccsdt.convergence", "convergence", reqs);
    this->addProduct("eomeeccsdt.R", "R", reqs);
}

template <typename U>
bool EOMEECCSDT<U>::run(TaskDAG& dag, const Arena& arena)
{
    auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const PointGroup& group = H.getABIJ().getGroup();
    int nirrep = group.getNumIrreps();

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->puttmp("D", new Denominator<U>(H));
    auto& TDAevecs = this->template get<vector<unique_vector<SpinorbitalTensor<U>>>>("TDAevecs");
    auto& TDAevals = this->template get<vector<vector<U>>>("TDAevals");

    this->puttmp("XMI", new SpinorbitalTensor<U>("X(mi)", H.getIJ()));
    this->puttmp("XAE", new SpinorbitalTensor<U>("X(ae)", H.getAB()));
    this->puttmp("XME", new SpinorbitalTensor<U>("X(me)", H.getIA()));
    this->puttmp("XABEJ", new SpinorbitalTensor<U>("X(ab,ej)", H.getABCI()));
    this->puttmp("XAMIJ", new SpinorbitalTensor<U>("X(am,ij)", H.getAIJK()));
    this->puttmp("XABEF", new SpinorbitalTensor<U>("X(ab,ef)", H.getABCD()));
    this->puttmp("XAMEI", new SpinorbitalTensor<U>("X(am,ei)", H.getAIBJ()));
    this->puttmp("XMNIJ", new SpinorbitalTensor<U>("X(mn,ij)", H.getIJKL()));

    vector<tuple<U,int,int,int>> tda_sorted;

    int tot_singlet = 0;
    int tot_triplet = 0;
    for (int i = 0;i < nirrep;i++)
    {
        vector<int> spin(TDAevals[i].size());

        for (int j = 0;j < TDAevals[i].size();j++)
        {
            if (scalar(TDAevecs[i][j]({1,0},{0,1})*TDAevecs[i][j]({0,0},{0,0})) < 0) {
                spin[j] = 1;
                tot_triplet++;
            }
            else
                tot_singlet++;
        }

        tda_sorted += zip(TDAevals[i],
                          spin,
                          vector<int>(TDAevals[i].size(), i),
                          vector<int>(range<int>(TDAevals[i].size())));
    }

    sort(tda_sorted);

    vector<vector<int>> root_idx(nirrep);
    if (ntriplet != 0)
    {
        ntriplet = min(ntriplet, tot_triplet);
        int nsinglet = min(tot_singlet, nroot-ntriplet);
        assert (nsinglet + ntriplet == nroot);
        int singlets_taken = 0;
        int triplets_taken = 0;
        int index = 0;

        while (singlets_taken + triplets_taken < nroot)
        {
            bool take_root = false;
            int this_spin = get<1>(tda_sorted[index]);
            if (this_spin == 0 and singlets_taken < nsinglet)
            {
                singlets_taken++;
                take_root = true;
            }
            if (this_spin == 1 && triplets_taken < ntriplet)
            {
                triplets_taken++;
                take_root = true;
            }
            if (take_root)
            {
                root_idx[get<2>(tda_sorted[index])].push_back(
                    get<3>(tda_sorted[index]));
                if (arena.rank == 0)
                    cout << "Took root " << index << " with spin "
                         << this_spin << endl;
            }
            index++;
        }
    }
    else
    {
        int roots_taken = 0;
        int index = 0;
        while (roots_taken < nroot)
        {
            int this_spin = get<1>(tda_sorted[index]);
            root_idx[get<2>(tda_sorted[index])].push_back(
                get<3>(tda_sorted[index]));
            roots_taken++;
            if (arena.rank == 0)
                cout << "Took root " << index << " with spin "
                     << this_spin << endl;
            index++;
        }
    }

    auto& Rs = this->puttmp("R", new unique_vector<ExcitationOperator<U,3>>());
    auto& Zs = this->puttmp("Z", new unique_vector<ExcitationOperator<U,3>>());

    for (int i = 0;i < nirrep;i++)
    {
        Rs.clear();
        Zs.clear();

        if (multiroot)
        {
            for (int j = 0;j < root_idx[i].size();j++)
            {
                Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(i));
                Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(i));
                ExcitationOperator<U,3>& R = Rs.back();
                R(0) = 0;
                R(1) = TDAevecs[i][root_idx[i][j]];
                R(2) = 0;
                R(3) = 0;
            }

            auto& davidson = this->puttmp("Davidson",
                new Davidson<ExcitationOperator<U,3>>(davidson_config, (int)root_idx[i].size()));

            previous.assign(root_idx[i].size(), numeric_limits<U>::max());

            Iterative<U>::run(dag, arena, root_idx[i].size());
        }
        else
        {
            auto& davidson = this->puttmp("Davidson",
                new Davidson<ExcitationOperator<U,3>>(davidson_config));

            for (int j = 0;j < root_idx[i].size();j++)
            {
                Logger::log(arena) << "Starting root number " << (j+1) << endl;

                Rs.clear();
                Zs.clear();

                Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(i));
                Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(i));
                ExcitationOperator<U,3>& R = Rs.back();
                R(0) = 0;
                R(1) = TDAevecs[i][root_idx[i][j]];
                R(2) = 0;
                R(3) = 0;

                R(1) = 0;
                vector<tkv_pair<U>> pairs = {{57, sqrt(0.5)}};
                R(1)({0,0},{0,0})({0,0}).writeRemoteData(pairs);
                R(1)({1,0},{0,1})({0,0}).writeRemoteData(pairs);

                triplet = scalar(R(1)({1,0},{0,1})*R(1)({0,0},{0,0})) < 0;

                if (triplet)
                {
                    Logger::log(arena) << "Triplet initial guess" << endl;
                }
                else
                {
                    Logger::log(arena)<< "Singlet initial guess" << endl;
                }

                previous.assign(1, numeric_limits<U>::max());

                Iterative<U>::run(dag, arena);

                if (!this->isConverged())
                {
                    this->error(arena) << "Root " << (j+1) << " did not converge." << endl;
                }

                davidson.getSolution(0, R);
                bool temp = scalar(R(1)({1,0},{0,1})*R(1)({0,0},{0,0})) < 0;
                if (temp)
                    this->log(arena) << "triplet solution found!" << endl;
                else
                    this->log(arena) << "singlet solution found!" << endl;
                if (triplet != temp)
                {
                    this->log(arena) << "WARNING: Spin character different from initial guess!" << endl;
                }

                davidson.nextRoot();
            }
        }
    }

    this->put("energy", new CTFTensor<U>("energy", arena, 1, {nroot}, {NS}, true));
    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void EOMEECCSDT<U>::iterate(const Arena& arena)
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

    auto& XMI   = this->template gettmp<SpinorbitalTensor<U>>("XMI"  );
    auto& XAE   = this->template gettmp<SpinorbitalTensor<U>>("XAE"  );
    auto& XME   = this->template gettmp<SpinorbitalTensor<U>>("XME"  );
    auto& XABEJ = this->template gettmp<SpinorbitalTensor<U>>("XABEJ");
    auto& XAMIJ = this->template gettmp<SpinorbitalTensor<U>>("XAMIJ");
    auto& XABEF = this->template gettmp<SpinorbitalTensor<U>>("XABEF");
    auto& XAMEI = this->template gettmp<SpinorbitalTensor<U>>("XAMEI");
    auto& XMNIJ = this->template gettmp<SpinorbitalTensor<U>>("XMNIJ");

    auto& D = this->template gettmp<Denominator<U>>("D");
    auto& davidson = this->template gettmp<Davidson<ExcitationOperator<U,3>>>("Davidson");

    auto& Rs = this->template gettmp<unique_vector<ExcitationOperator<U,3>>>("R");
    auto& Zs = this->template gettmp<unique_vector<ExcitationOperator<U,3>>>("Z");

    for (int root = 0;root < this->nsolution();root++)
    {
        ExcitationOperator<U,3>& R = Rs[root];
        ExcitationOperator<U,3>& Z = Zs[root];
        Z = 0;

        if (triplet)
        {
            0.5*R(1)({0,0},{0,0})[  "ai"] -= 0.5*R(1)({1,0},{0,1})[  "ai"];
                R(1)({1,0},{0,1})[  "ai"]  =    -R(1)({0,0},{0,0})[  "ai"];
            0.5*R(2)({1,0},{0,1})["abij"] -= 0.5*R(2)({1,0},{0,1})["baji"];
            0.5*R(2)({0,0},{0,0})["abij"] -= 0.5*R(2)({2,0},{0,2})["abij"];
                R(2)({2,0},{0,2})["abij"]  =    -R(2)({0,0},{0,0})["abij"];
        }
        else
        {
            0.5*R(1)({0,0},{0,0})[  "ai"] += 0.5*R(1)({1,0},{0,1})[  "ai"];
                R(1)({1,0},{0,1})[  "ai"]  =     R(1)({0,0},{0,0})[  "ai"];
            0.5*R(2)({1,0},{0,1})["abij"] += 0.5*R(2)({1,0},{0,1})["baji"];
            0.5*R(2)({0,0},{0,0})["abij"] += 0.5*R(2)({2,0},{0,2})["abij"];
                R(2)({2,0},{0,2})["abij"]  =     R(2)({0,0},{0,0})["abij"];
        }

          XME[    "me"]  =      WMNEF["mnef"]*R(1)[    "fn"];

          XMI[    "mi"]  =      WMNEJ["nmei"]*R(1)[    "en"];
          XMI[    "mi"] +=  0.5*WMNEF["mnef"]*R(2)[  "efin"];

          XAE[    "ae"]  =      WAMEF["amef"]*R(1)[    "fm"];
          XAE[    "ae"] -=  0.5*WMNEF["mnef"]*R(2)[  "afmn"];

        XABEF[  "abef"]  =  0.5*WMNEF["mnef"]*R(2)[  "abmn"];
        XABEF[  "abef"] -=      WAMEF["amef"]*R(1)[    "bm"];

        XAMEI[  "amei"]  =      WMNEF["mnef"]*R(2)[  "afni"];
        XAMEI[  "amei"] +=      WAMEF["amef"]*R(1)[    "fi"];
        XAMEI[  "amei"] -=      WMNEJ["nmei"]*R(1)[    "an"];

        XMNIJ[  "mnij"]  =  0.5*WMNEF["mnef"]*R(2)[  "efij"];
        XMNIJ[  "mnij"] +=      WMNEJ["mnej"]*R(1)[    "ei"];

        XAMIJ[  "amij"]  =  0.5*WAMEF["amef"]*R(2)[  "efij"];
        XAMIJ[  "amij"] +=      WMNEJ["nmej"]*R(2)[  "aein"];
        XAMIJ[  "amij"] +=      WAMEI["amej"]*R(1)[    "ei"];
        XAMIJ[  "amij"] -=      WMNIJ["nmij"]*R(1)[    "an"];
        XAMIJ[  "amij"] +=        XME[  "me"]*T(2)[  "aeij"];
        XAMIJ[  "amij"] +=  0.5*WMNEF["mnef"]*R(3)["aefijn"];

        XABEJ[  "abej"]  =  0.5*WMNEJ["mnej"]*R(2)[  "abmn"];
        XABEJ[  "abej"] +=      WAMEF["amef"]*R(2)[  "fbmj"];
        XABEJ[  "abej"] +=      WABEF["abef"]*R(1)[    "fj"];
        XABEJ[  "abej"] -=      WAMEI["amej"]*R(1)[    "bm"];
        XABEJ[  "abej"] -=  0.5*WMNEF["mnef"]*R(3)["afbmnj"];

         Z(1)[    "ai"]  =        FAE[  "ae"]*R(1)[    "ei"];
         Z(1)[    "ai"] -=        FMI[  "mi"]*R(1)[    "am"];
         Z(1)[    "ai"] -=      WAMEI["amei"]*R(1)[    "em"];
         Z(1)[    "ai"] +=        FME[  "me"]*R(2)[  "aeim"];
         Z(1)[    "ai"] +=  0.5*WAMEF["amef"]*R(2)[  "efim"];
         Z(1)[    "ai"] -=  0.5*WMNEJ["mnei"]*R(2)[  "eamn"];
         Z(1)[    "ai"] += 0.25*WMNEF["mnef"]*R(3)["efamni"];

         Z(2)[  "abij"]  =      WABEJ["abej"]*R(1)[    "ei"];
         Z(2)[  "abij"] -=      WAMIJ["amij"]*R(1)[    "bm"];
         Z(2)[  "abij"] +=        FAE[  "ae"]*R(2)[  "ebij"];
         Z(2)[  "abij"] -=        FMI[  "mi"]*R(2)[  "abmj"];
         Z(2)[  "abij"] +=        XAE[  "ae"]*T(2)[  "ebij"];
         Z(2)[  "abij"] -=        XMI[  "mi"]*T(2)[  "abmj"];
         Z(2)[  "abij"] +=  0.5*WMNIJ["mnij"]*R(2)[  "abmn"];
         Z(2)[  "abij"] +=  0.5*WABEF["abef"]*R(2)[  "efij"];
         Z(2)[  "abij"] -=      WAMEI["amei"]*R(2)[  "ebmj"];
         Z(2)[  "abij"] +=        FME[  "me"]*R(3)["abeijm"];
         Z(2)[  "abij"] +=        XME[  "me"]*T(3)["abeijm"];
         Z(2)[  "abij"] +=  0.5*WAMEF["amef"]*R(3)["efbimj"];
         Z(2)[  "abij"] -=  0.5*WMNEJ["mnej"]*R(3)["aebimn"];

          XMI[    "mi"] +=        FME[  "me"]*R(1)[    "ei"];

          XAE[    "ae"] -=        FME[  "me"]*R(1)[    "am"];

         Z(3)["abcijk"]  =      WABEJ["abej"]*R(2)[  "ecik"];
         Z(3)["abcijk"] -=      WAMIJ["amij"]*R(2)[  "bcmk"];
         Z(3)["abcijk"] +=      XABEJ["abej"]*T(2)[  "ecik"];
         Z(3)["abcijk"] -=      XAMIJ["amij"]*T(2)[  "bcmk"];
         Z(3)["abcijk"] +=        FAE[  "ae"]*R(3)["ebcijk"];
         Z(3)["abcijk"] -=        FMI[  "mi"]*R(3)["abcmjk"];
         Z(3)["abcijk"] +=        XAE[  "ae"]*T(3)["ebcijk"];
         Z(3)["abcijk"] -=        XMI[  "mi"]*T(3)["abcmjk"];
         Z(3)["abcijk"] +=  0.5*WABEF["abef"]*R(3)["efcijk"];
         Z(3)["abcijk"] -=      WAMEI["amei"]*R(3)["ebcmjk"];
         Z(3)["abcijk"] +=  0.5*WMNIJ["mnij"]*R(3)["abcmnk"];
         Z(3)["abcijk"] +=  0.5*XABEF["abef"]*T(3)["efcijk"];
         Z(3)["abcijk"] -=      XAMEI["amei"]*T(3)["ebcmjk"];
         Z(3)["abcijk"] +=  0.5*XMNIJ["mnij"]*T(3)["abcmnk"];

         printf("R3: %.15f\n", scalar(R(3)*R(3))/36);
         printf("R1: %.15f\n", scalar(R(1)*R(1)));
         printf("R2: %.15f\n", scalar(R(2)*R(2))/4);

         printf("XME:   %.15f\n", XME  ({0,0},{0,0}).norm(2));
         printf("XAE:   %.15f\n", XAE  ({0,0},{0,0}).norm(2));
         printf("XMI:   %.15f\n", XMI  ({0,0},{0,0}).norm(2));
         printf("XABEJ: %.15f\n", XABEJ({1,0},{1,0}).norm(2));
         printf("XAMIJ: %.15f\n", XAMIJ({0,1},{0,1}).norm(2));
         printf("XABEF: %.15f\n", XABEF({1,0},{1,0}).norm(2));
         printf("XMAEI: %.15f\n", XAMEI({1,0},{0,1}).norm(2));
         printf("XAMEI: %.15f\n", XAMEI({0,1},{0,1}).norm(2));
         printf("XMNIJ: %.15f\n", XMNIJ({0,1},{0,1}).norm(2));

         printf("FME:   %.15f\n", FME  ({0,0},{0,0}).norm(2));
         printf("FAE:   %.15f\n", FAE  ({0,0},{0,0}).norm(2));
         printf("FMI:   %.15f\n", FMI  ({0,0},{0,0}).norm(2));
         printf("WABEJ: %.15f\n", WABEJ({1,0},{1,0}).norm(2));
         printf("WAMIJ: %.15f\n", WAMIJ({0,1},{0,1}).norm(2));
         printf("WABEF: %.15f\n", WABEF({1,0},{1,0}).norm(2));
         printf("WMAEI: %.15f\n", WAMEI({1,0},{0,1}).norm(2));
         printf("WAMEI: %.15f\n", WAMEI({0,1},{0,1}).norm(2));
         printf("WMNIJ: %.15f\n", WMNIJ({0,1},{0,1}).norm(2));

         printf("Z3: %.15f\n", scalar(Z(3)*Z(3))/36);
         printf("Z1: %.15f\n", scalar(Z(1)*Z(1)));
         printf("Z2: %.15f\n", scalar(Z(2)*Z(2))/4);
    }

    vector<U> energies = davidson.extrapolate(Rs, Zs, D);

    for (int i = 0;i < this->nsolution();i++)
    {
        this->energy(i) = energies[i];
        this->conv(i) = max(aquarius::abs(energies[i]-previous[i]), Zs[i].norm(00));
        previous[i] = energies[i];
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
            int 10,
    jacobi?
            bool false,
    num_reduce?
            int 3,
    compaction?
            enum { discrete, continuous },
}

)";

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSDT);
REGISTER_TASK(aquarius::cc::EOMEECCSDT<double>, "eomeeccsdt",spec);
