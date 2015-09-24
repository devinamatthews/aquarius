#include "eomeeccsd.hpp"

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
EOMEECCSD<U>::EOMEECCSD(const string& name, Config& config)
: Iterative<U>(name, config),
  davidson_config(config.get("davidson")),
  nroot(config.get<int>("nroot")),
  ntriplet(config.get<int>("ntriplet")),
  nsinglet(config.get<int>("nsinglet")),
  multiroot(config.get<bool>("multiroot"))
{
    vector<Requirement> reqs;
    reqs.emplace_back("ccsd.T", "T");
    reqs.emplace_back("ccsd.Hbar", "Hbar");
    reqs.emplace_back("tda.TDAevals", "TDAevals");
    reqs.emplace_back("tda.TDAevecs", "TDAevecs");
    this->addProduct("eomeeccsd.energy", "energy", reqs);
    this->addProduct("eomeeccsd.convergence", "convergence", reqs);
    this->addProduct("eomeeccsd.R", "R", reqs);
}

template <typename U>
bool EOMEECCSD<U>::run(TaskDAG& dag, const Arena& arena)
{
    auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

    const PointGroup& group = H.getABIJ().getGroup();
    int nirrep = group.getNumIrreps();

    const Space& occ = H.occ;
    const Space& vrt = H.vrt;

    this->puttmp("D", new Denominator<U>(H));
    auto& TDAevecs = this->template get<vector<unique_vector<SpinorbitalTensor<U>>>>("TDAevecs");
    auto& TDAevals = this->template get<vector<vector<U>>>("TDAevals");

    vector<tuple<U,int,int,int>> tda_sorted;

    int tot_singlet = 0;
    int tot_triplet = 0;
    for (int i = 0;i < nirrep;i++)
    {
        vector<int> spin(TDAevals[i].size());

        for (int j = 0;j < TDAevals[i].size();j++)
        {
            auto& R1a = TDAevecs[i][j]({0,0},{0,0});
            auto& R1b = TDAevecs[i][j]({1,0},{0,1});
            U olap = scalar(R1a*R1b);
            //cout << i << " " << j << " " << setprecision(15) << olap << " " << TDAevals[i][j] << endl;
            if (olap < 0)
            {
                spin[j] = 1;
                R1a -=  R1b;
                R1b  = -R1a;
                tot_triplet++;
            }
            else
            {
                R1a += R1b;
                R1b  = R1a;
                tot_singlet++;
            }
            U nrm = sqrt(2*scalar(R1a*R1a));
            R1a /= nrm;
            R1b /= nrm;
        }

        tda_sorted += zip(TDAevals[i],
                          spin,
                          vector<int>(TDAevals[i].size(), i),
                          vector<int>(range<int>(TDAevals[i].size())));
    }

    assert(nroot > 0);
    assert(nsinglet >= 0);
    assert(ntriplet >= 0);
    assert(nsinglet+ntriplet <= nroot);
    assert(nsinglet <= tot_singlet);
    assert(ntriplet <= tot_triplet);
    assert(nroot <= tot_singlet+tot_triplet);

    sort(tda_sorted);

    vector<tuple<U,int,int,int>> roots;

    for (auto& root : tda_sorted)
    {
        if (nsinglet == 0) break;
        if (get<1>(root) == 0)
        {
            nsinglet--;
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    for (auto& root : tda_sorted)
    {
        if (ntriplet == 0) break;
        if (get<1>(root) == 1)
        {
            ntriplet--;
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    for (auto& root : tda_sorted)
    {
        if (nroot == 0) break;
        if (get<1>(root) != -1)
        {
            nroot--;
            roots.push_back(root);
            get<1>(root) = -1;
        }
    }

    auto& Rs = this->puttmp("R", new unique_vector<ExcitationOperator<U,2>>());
    //auto& Vs = this->puttmp("V", new unique_vector<ExcitationOperator<U,2>>());
    auto& Zs = this->puttmp("Z", new unique_vector<ExcitationOperator<U,2>>());

    int idx = 1;
    for (int irrep = 0;irrep < nirrep;irrep++)
    {
        for (int spin : {0,1})
        {
            triplet = spin == 1;

            //Vs.clear();
            Rs.clear();
            Zs.clear();

            if (multiroot)
            {
                /*
                for (int j = 0;j < root_idx[i].size();j++)
                {
                    Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(i));
                    Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(i));
                    ExcitationOperator<U,2>& R = Rs.back();
                    R(0) = 0;
                    R(1) = TDAevecs[i][root_idx[i][j]];
                    R(2) = 0;
                }

                auto& davidson = this->puttmp("Davidson",
                    new Davidson<ExcitationOperator<U,2>>(davidson_config, (int)root_idx[i].size()));

                previous.assign(root_idx[i].size(), numeric_limits<U>::max());

                Iterative<U>::run(dag, arena, root_idx[i].size());
                */

                //TODO: cleanup and retest
                assert(0);

                /*
                for (int j = 0;j < root_idx[i].size();j++)
                {
                    if (this->isConverged(j))
                    {
                        Vs.emplace_back("V", arena, occ, vrt, group.getIrrep(i));
                        ExcitationOperator<U,2>& V = Vs.back();
                        davidson.getSolution(j, V);
                        V /= sqrt(aquarius::abs(scalar(conj(V)*V)));
                    }
                }
                */
            }
            else
            {
                auto& davidson = this->puttmp("Davidson",
                    new Davidson<ExcitationOperator<U,2>>(davidson_config));

                for (auto& root : roots)
                {
                    if (spin != get<1>(root)) continue;
                    if (irrep != get<2>(root)) continue;

                    int which = get<3>(root);

                    Logger::log(arena) << "Starting root number " << idx << endl;
                    Logger::log(arena) << "Guess energy: " << fixed << setprecision(12) << get<0>(root) << endl;

                    Rs.clear();
                    Zs.clear();

                    Rs.emplace_back("R", arena, occ, vrt, group.getIrrep(irrep));
                    Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(irrep));

                    this->puttmp("XMI", new SpinorbitalTensor<U>("X(mi)", arena, group, group.getIrrep(irrep), {vrt,occ}, {0,1}, {0,1}));
                    this->puttmp("XAE", new SpinorbitalTensor<U>("X(ae)", arena, group, group.getIrrep(irrep), {vrt,occ}, {1,0}, {1,0}));

                    ExcitationOperator<U,2>& R = Rs.back();
                    R(0) = 0;
                    R(1) = TDAevecs[irrep][which];
                    R(2) = 0;

                    if (triplet)
                    {
                        Logger::log(arena) << "Triplet initial guess" << endl;
                    }
                    else
                    {
                        Logger::log(arena)<< "Singlet initial guess" << endl;
                    }

                    bool print_vecs;
                    print_vecs = false;

                    if (print_vecs)
                    {
                        vector<U> temp1;
                        vector<U> temp2;
                        R(1)({1,0},{0,1})({0,0}).getAllData(temp1);
                        R(1)({0,0},{0,0})({0,0}).getAllData(temp2);

                        if (arena.rank == 0)
                        {
                            cout << " " << endl;
                            cout << "Root " << idx << " R1" << endl;
                            for (int ii=0; ii<temp1.size(); ii++)
                            {
                                cout << ii << " " << temp1[ii] << " " << temp2[ii] << endl;
                            }
                            cout << " " << endl;
                        }
                    }

                    previous.assign(1, numeric_limits<U>::max());

                    Iterative<U>::run(dag, arena);

                    if (!this->isConverged())
                    {
                        this->error(arena) << "Root " << idx << " did not converge." << endl;
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

                    if (print_vecs)
                    {
                        vector<U> temp1;
                        vector<U> temp2;
                        R(1)({1,0},{0,1})({0,0}).getAllData(temp1);
                        R(1)({0,0},{0,0})({0,0}).getAllData(temp2);
                        vector<tuple<U,U,U,int>> amps_sorted;

                        for (int ii=0; ii < temp1.size(); ii++)
                            amps_sorted.emplace_back(-abs(temp1[ii]),temp1[ii],temp2[ii],ii);

                        sort(amps_sorted);

                        double norm = sqrt(aquarius::abs(scalar(conj(R)*R)));

                        if (arena.rank == 0)
                        {
                            cout << " " << endl;
                            cout << "sqrt(aquarius::abs(scalar(conj(R)*R))) = " << norm << endl;
                            cout << " " << endl;
                            cout << "Root " << idx << " R1" << endl;
                            for (int ii=0; ii<30; ii++)
                            {
                                cout << get<3>(amps_sorted[ii]) << " " << get<1>(amps_sorted[ii]) << " " << get<2>(amps_sorted[ii]) << endl;
                            }
                            cout << " " << endl;
                        }
                    }

                    davidson.nextRoot();
                    idx++;
                }
            }
        }
    }

    this->put("energy", new CTFTensor<U>("energy", arena, 1, {nroot}, {NS}, true));
    this->put("convergence", new U(this->conv()));

    return true;
}

template <typename U>
void EOMEECCSD<U>::iterate(const Arena& arena)
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
    //auto& Vs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("V");
    auto& Zs = this->template gettmp<unique_vector<ExcitationOperator<U,2>>>("Z");

    for (int root = 0;root < this->nsolution();root++)
    {
        ExcitationOperator<U,2>& R = Rs[root];
        ExcitationOperator<U,2>& Z = Zs[root];

        double sign = triplet ? -1 : 1;

        0.5*R(1)({0,0},{0,0})[  "ai"] += 0.5*sign*R(1)({1,0},{0,1})[  "ai"];
            R(1)({1,0},{0,1})[  "ai"]  =     sign*R(1)({0,0},{0,0})[  "ai"];

        0.5*R(2)({1,0},{0,1})["abij"] += 0.5*sign*R(2)({1,0},{0,1})["baji"];
        0.5*R(2)({0,0},{0,0})["abij"] += 0.5*sign*R(2)({2,0},{0,2})["abij"];
            R(2)({2,0},{0,2})["abij"]  =     sign*R(2)({0,0},{0,0})["abij"];

         XMI[  "mi"]  =     WMNEJ["nmei"]*R(1)[  "en"];
         XMI[  "mi"] += 0.5*WMNEF["mnef"]*R(2)["efin"];
         XAE[  "ae"]  =     WAMEF["amef"]*R(1)[  "fm"];
         XAE[  "ae"] -= 0.5*WMNEF["mnef"]*R(2)["afmn"];

        Z(1)[  "ai"]  =       FAE[  "ae"]*R(1)[  "ei"];
        Z(1)[  "ai"] -=       FMI[  "mi"]*R(1)[  "am"];
        Z(1)[  "ai"] -=     WAMEI["amei"]*R(1)[  "em"];
        Z(1)[  "ai"] +=       FME[  "me"]*R(2)["aeim"];
        Z(1)[  "ai"] += 0.5*WAMEF["amef"]*R(2)["efim"];
        Z(1)[  "ai"] -= 0.5*WMNEJ["mnei"]*R(2)["eamn"];

        Z(2)["abij"]  =     WABEJ["abej"]*R(1)[  "ei"];
        Z(2)["abij"] -=     WAMIJ["amij"]*R(1)[  "bm"];
        Z(2)["abij"] +=       FAE[  "ae"]*R(2)["ebij"];
        Z(2)["abij"] -=       FMI[  "mi"]*R(2)["abmj"];
        Z(2)["abij"] +=       XAE[  "ae"]*T(2)["ebij"];
        Z(2)["abij"] -=       XMI[  "mi"]*T(2)["abmj"];
        Z(2)["abij"] += 0.5*WMNIJ["mnij"]*R(2)["abmn"];
        Z(2)["abij"] += 0.5*WABEF["abef"]*R(2)["efij"];
        Z(2)["abij"] -=     WAMEI["amei"]*R(2)["ebmj"];
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
    int 1,
nsinglet?
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::EOMEECCSD);
REGISTER_TASK(aquarius::cc::EOMEECCSD<double>, "eomeeccsd",spec);
