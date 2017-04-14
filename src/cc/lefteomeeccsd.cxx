#include "util/global.hpp"

#include "convergence/davidson.hpp"
#include "util/iterative.hpp"
#include "operator/2eoperator.hpp"
#include "operator/st2eoperator.hpp"
#include "operator/excitationoperator.hpp"
#include "operator/denominator.hpp"

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
class LeftEOMEECCSD : public Iterative<U>
{
    protected:
        Config davidson_config;
        bool multiroot;
        bool triplet = false;
        vector<U> previous;

    public:
        LeftEOMEECCSD(const string& name, Config& config)
        : Iterative<U>(name, config),
          davidson_config(config.get("davidson")),
          multiroot(config.get<bool>("multiroot"))
        {
            vector<Requirement> reqs;
            reqs.emplace_back("ccsd.T", "T");
            reqs.emplace_back("ccsd.Hbar", "Hbar");
            reqs.emplace_back("eomeeccsd.energy", "energy");
            reqs.emplace_back("eomeeccsd.R", "R");
            this->addProduct("eomeeccsd.L", "L", reqs);
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            auto& H = this->template get<STTwoElectronOperator<U>>("Hbar");

            const PointGroup& group = H.getABIJ().getGroup();
            int nirrep = group.getNumIrreps();

            const Space& occ = H.occ;
            const Space& vrt = H.vrt;

            this->puttmp("D", new Denominator<U>(H));

            auto& Es = this->template get<vector<U>>("energy");
            auto& Rs = this->template get<unique_vector<ExcitationOperator<U,2>>>("R");

            auto& Ls = this->put   ("L", new unique_vector<DeexcitationOperator<U,2>>());
            auto& Cs = this->puttmp("C", new unique_vector<DeexcitationOperator<U,2>>());
            auto& Zs = this->puttmp("Z", new unique_vector<DeexcitationOperator<U,2>>());

            int idx = 1;
            for (int irrep = 0;irrep < nirrep;irrep++)
            {
                Cs.clear();
                Zs.clear();

                if (multiroot)
                {
                    assert(0);
                }
                else
                {
                    auto& davidson = this->puttmp("Davidson",
                        new Davidson<DeexcitationOperator<U,2>>(davidson_config));

                    for (int root = 0;root < Rs.size();root++)
                    {
                        auto& R = Rs[root];
                        U E = Es[root];

                        if (!(R.getRepresentation()*group.getIrrep(irrep))) continue;

                        davidson.reset(E);

                        Logger::log(arena) << "Starting root number " << idx << endl;
                        Logger::log(arena) << "Guess energy: " << fixed << setprecision(12) << E << endl;

                        Cs.clear();
                        Zs.clear();

                        Cs.emplace_back("C", arena, occ, vrt, group.getIrrep(irrep));
                        Zs.emplace_back("Z", arena, occ, vrt, group.getIrrep(irrep));

                        this->puttmp("XIM", new SpinorbitalTensor<U>("X(im)", arena, group, group.getIrrep(irrep), {vrt,occ}, {0,1}, {0,1}));
                        this->puttmp("XEA", new SpinorbitalTensor<U>("X(ea)", arena, group, group.getIrrep(irrep), {vrt,occ}, {1,0}, {1,0}));

                        DeexcitationOperator<U,2>& L = Cs.back();
                        L(0) = 0;
                        L(1)["ia"] = R(1)["ai"];
                        L(2)["ijab"] = R(2)["abij"];
                        L /= sqrt(aquarius::abs(scalar(conj(L)*L)));

                        triplet = scalar(L(1)({0,1},{1,0})*L(1)({0,0},{0,0})) < 0;

                        bool print_vecs;
                        print_vecs = false;

                        if (print_vecs)
                        {
                            vector<U> temp1;
                            vector<U> temp2;
                            L(1)({0,1},{1,0})({0,0}).getAllData(temp1);
                            L(1)({0,0},{0,0})({0,0}).getAllData(temp2);

                            if (arena.rank == 0)
                            {
                                cout << " " << endl;
                                cout << "Root " << idx << " L1" << endl;
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

                        davidson.getSolution(0, L);
                        bool temp = scalar(L(1)({0,1},{1,0})*L(1)({0,0},{0,0})) < 0;
                        if (temp)
                            this->log(arena) << "triplet solution found!" << endl;
                        else
                            this->log(arena) << "singlet solution found!" << endl;
                        if (triplet != temp)
                        {
                            this->log(arena) << "WARNING: Spin character different from initial guess!" << endl;
                        }

                        L /= scalar(conj(L(1)["ia"])*R(1)["ai"]) +
                             0.25*scalar(conj(L(2)["ijab"])*R(2)["abij"]);

                        Ls.emplace_back(L);

                        if (print_vecs)
                        {
                            vector<U> temp1;
                            vector<U> temp2;
                            L(1)({0,1},{1,0})({0,0}).getAllData(temp1);
                            L(1)({0,0},{0,0})({0,0}).getAllData(temp2);
                            vector<tuple<U,U,U,int>> amps_sorted;

                            for (int ii=0; ii < temp1.size(); ii++)
                                amps_sorted.emplace_back(-abs(temp1[ii]),temp1[ii],temp2[ii],ii);

                            sort(amps_sorted);

                            double norm = sqrt(aquarius::abs(scalar(conj(L)*L)));

                            if (arena.rank == 0)
                            {
                                cout << " " << endl;
                                cout << "sqrt(aquarius::abs(scalar(conj(L)*L))) = " << norm << endl;
                                cout << " " << endl;
                                cout << "Root " << idx << " L1" << endl;
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

            auto& T = this->template get<ExcitationOperator<U,2>>("T");

            auto& XIM = this->template gettmp<SpinorbitalTensor<U>>("XIM");
            auto& XEA = this->template gettmp<SpinorbitalTensor<U>>("XEA");

            auto& D = this->template gettmp<Denominator<U>>("D");
            auto& davidson = this->template gettmp<Davidson<DeexcitationOperator<U,2>>>("Davidson");

            auto& Cs = this->template gettmp<unique_vector<DeexcitationOperator<U,2>>>("C");
            auto& Zs = this->template gettmp<unique_vector<DeexcitationOperator<U,2>>>("Z");

            for (int root = 0;root < this->nsolution();root++)
            {
                DeexcitationOperator<U,2>& L = Cs[root];
                DeexcitationOperator<U,2>& Z = Zs[root];

                double sign = triplet ? -1 : 1;

                0.5*L(1)({0,0},{0,0})[  "ia"] += 0.5*sign*L(1)({0,1},{1,0})[  "ia"];
                    L(1)({0,1},{1,0})[  "ia"]  =     sign*L(1)({0,0},{0,0})[  "ia"];

                0.5*L(2)({0,1},{1,0})["ijab"] += 0.5*sign*L(2)({0,1},{1,0})["jiba"];
                0.5*L(2)({0,0},{0,0})["ijab"] += 0.5*sign*L(2)({0,2},{2,0})["ijab"];
                    L(2)({0,2},{2,0})["ijab"]  =     sign*L(2)({0,0},{0,0})["ijab"];

                 XIM[  "im"]  =  0.5*T(2)["efmn"]* L(2)["inef"];
                 XEA[  "ea"]  = -0.5*T(2)["efmn"]* L(2)["mnaf"];

                Z(1)[  "ia"]  =       FAE[  "ea"]* L(1)[  "ie"];
                Z(1)[  "ia"] -=       FMI[  "im"]* L(1)[  "ma"];
                Z(1)[  "ia"] -=     WAMEI["eiam"]* L(1)[  "me"];
                Z(1)[  "ia"] -=       XEA[  "ef"]*WAMEF["fiea"];
                Z(1)[  "ia"] -=       XIM[  "nm"]*WMNEJ["iman"];
                Z(1)[  "ia"] += 0.5*WABEJ["efam"]* L(2)["imef"];
                Z(1)[  "ia"] -= 0.5*WAMIJ["eimn"]* L(2)["mnea"];

                Z(2)["ijab"]  =       FME[  "ia"]* L(1)[  "jb"];
                Z(2)["ijab"] +=     WAMEF["ejab"]* L(1)[  "ie"];
                Z(2)["ijab"] -=     WMNEJ["ijam"]* L(1)[  "mb"];
                Z(2)["ijab"] +=       FAE[  "ea"]* L(2)["ijeb"];
                Z(2)["ijab"] -=       FMI[  "im"]* L(2)["mjab"];
                Z(2)["ijab"] +=       XEA[  "ea"]*WMNEF["ijeb"];
                Z(2)["ijab"] -=       XIM[  "im"]*WMNEF["mjab"];
                Z(2)["ijab"] += 0.5*WMNIJ["ijmn"]* L(2)["mnab"];
                Z(2)["ijab"] += 0.5*WABEF["efab"]* L(2)["ijef"];
                Z(2)["ijab"] -=     WAMEI["eiam"]* L(2)["mjeb"];
            }

            vector<U> energies = davidson.extrapolate(Cs, Zs, D);

            for (int i = 0;i < this->nsolution();i++)
            {
                this->energy(i) = energies[i];
                this->conv(i) = Zs[i].norm(00);
            }
        }
};

}
}

static const char* spec = R"(

multiroot?
    bool false,
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

INSTANTIATE_SPECIALIZATIONS(aquarius::cc::LeftEOMEECCSD);
REGISTER_TASK(aquarius::cc::LeftEOMEECCSD<double>, "lefteomeeccsd",spec);
