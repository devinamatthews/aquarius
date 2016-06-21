#ifndef _AQUARIUS_TASKS_JELLIUM_HPP_
#define _AQUARIUS_TASKS_JELLIUM_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/tensor.hpp"
#include "frameworks/logging.hpp"
#include "frameworks/operator.hpp"

using namespace aquarius::task;
using namespace aquarius::logging;
using namespace aquarius::symmetry;
using namespace aquarius::tensor;
using namespace aquarius::op;

namespace aquarius
{
namespace jellium
{

class Jellium : public task::Task
{
    protected:
        int nelec;
        int norb;
        double radius;
        vector<vec3> gvecs;
        int nocc;
        double V;
        double L;
        double PotVm;

        void writeIntegrals(bool pvirt, bool qvirt, bool rvirt, bool svirt,
                            Tensor<SPINORBITAL> tensor)
        {
            vector<vector<int>> all_spins;
            if (pvirt == qvirt && rvirt == svirt)
                all_spins = {{1,1,1,1}, {1,0,1,0}, {0,0,0,0}};
            else if (pvirt == qvirt)
                all_spins = {{1,1,1,1}, {1,0,1,0}, {1,0,0,1}, {0,0,0,0}};
            else if (rvirt == svirt)
                all_spins = {{1,1,1,1}, {1,0,1,0}, {0,1,1,0}, {0,0,0,0}};
            else
                all_spins = {{1,1,1,1}, {1,0,1,0}, {1,0,0,1}, {0,1,1,0}, {0,1,0,1}, {0,0,0,0}};

            for (auto& spins : all_spins)
            {
                KeyValueVector pairs;
                tensor.getLocalDataBySpin(spins, pairs);

                int np = (pvirt ? norb-nocc : nocc);
                int nq = (qvirt ? norb-nocc : nocc);
                int nr = (rvirt ? norb-nocc : nocc);
                int ns = (svirt ? norb-nocc : nocc);

                for (size_t i = 0;i < pairs.size();i++)
                {
                    auto k = pairs.key(i);
                    int p = k%np;
                    k /= np;
                    int q = k%nq;
                    k /= nq;
                    int r = k%nr;
                    k /= nr;
                    int s = k;

                    if (pvirt) p += nocc;
                    if (qvirt) q += nocc;
                    if (rvirt) r += nocc;
                    if (svirt) s += nocc;

                    vec3 pr = gvecs[p]-gvecs[r];
                    vec3 ps = gvecs[p]-gvecs[s];
                    vec3 qr = gvecs[q]-gvecs[r];
                    vec3 qs = gvecs[q]-gvecs[s];

                    double val = 0;

                    if (spins[0] == spins[2] && norm2(pr+qs) < 1e-12)
                    {
                        val += (p == r ? PotVm : 1/(M_PI*L*norm2(pr)));
                    }
                    if (spins[0] == spins[3] && norm2(ps+qr) < 1e-12)
                    {
                        val -= (p == s ? PotVm : 1/(M_PI*L*norm2(ps)));
                    }

                    pairs.value(i, val);
                }

                tensor.setLocalDataBySpin(spins, pairs);
            }
        }

    public:
        Jellium(const string& name, Config& config)
        : Task(name, config),
          nelec(config.get<int>("num_electrons")),
          norb(config.get<int>("num_orbitals")),
          radius(config.get<double>("radius"))
        {
            vector<Requirement> reqs;
            addProduct("scf.energy", "energy", reqs);
            addProduct("scf.E", "E", reqs);
            addProduct("scf.F", "F", reqs);
            addProduct("scf.D", "D", reqs);
            addProduct("moints", "H", reqs);

            int d = config.get<int>("dimension");
            assert(d == 3);
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            vector<double> glen;

            for (int r = 0;;r++)
            {
                int n = 0;
                glen.clear();
                gvecs.clear();
                for (int x = -r;x <= r;x++)
                for (int y = -r;y <= r;y++)
                for (int z = -r;z <= r;z++)
                {
                    if (sqrt(x*x+y*y+z*z) < r)
                    {
                        glen.push_back(sqrt(x*x+y*y+z*z));
                        gvecs.emplace_back(x, y, z);
                        n++;
                    }
                }

                if (n >= norb) break;
            }

            cosort(glen, gvecs);

            assert(norb == gvecs.size() || std::abs(glen[norb-1] - glen[norb]) > 1e-12);
            gvecs.resize(norb);

            nocc = nelec/2;
            assert(0 < nocc && nocc <= norb);
            assert(nocc == norb || std::abs(glen[nocc-1] - glen[nocc]) > 1e-12);

            V = nelec*(4.0/3.0)*M_PI*pow(radius,3);
            L = pow(V, 1.0/3.0);
            PotVm = 2.83729747948149/L;

            vector<vector<double>> E = {vector<double>(norb)};

            int nvrt = norb-nocc;

            for (int i = 0;i < norb;i++)
            {
                E[0][i] = 2*(M_PI/L)*(M_PI/L)*norm2(gvecs[i]);
                for (int j = 0;j < nocc;j++)
                {
                    if (i == j)
                    {
                        E[0][i] -= PotVm;
                    }
                    else
                    {
                        E[0][i] -= 1/(M_PI*L*norm2(gvecs[i]-gvecs[j]));
                    }
                }
            }

            put("E", vector<vector<vector<double>>>{E,E});

            double energy = 0;
            for (int i = 0;i < nocc;i++)
            {
                energy += 2*E[0][i];
                for (int j = 0;j < nocc;j++)
                {
                    if (i == j)
                    {
                        energy += PotVm;
                    }
                    else
                    {
                        energy += 1/(M_PI*L*norm2(gvecs[i]-gvecs[j]));
                    }
                }
            }

            Tensor<SPINORBITAL> F = put<Tensor<>>("F",
                Tensor<SPINORBITAL>::construct("F", SOInit({norb}, {norb}, {1}, {1})));
            Tensor<SPINORBITAL> D = put<Tensor<>>("D", F.construct("D"));
            put("energy", energy);

            Logger::log(arena) << "SCF energy = " << setprecision(15) << energy << endl;

            KeyValueVector dpairs, fpairs;
            for (int i = 0;i < nocc;i++)
            {
                dpairs.push_back(i*norb+i, 1.0);
            }
            for (int i = 0;i < norb;i++)
            {
                fpairs.push_back(i*norb+i, E[0][i]);
            }

            D.setDataBySpin({1,1}, dpairs);
            D.setDataBySpin({0,0}, dpairs);
            F.setDataBySpin({1,1}, fpairs);
            F.setDataBySpin({0,0}, fpairs);

            auto& H = put("H", new TwoElectronOperator("H", SOInit({nvrt,nocc}, {nvrt,nocc})));

            KeyValueVector abpairs, ijpairs;
            for (int i = 0;i < nocc;i++)
            {
                ijpairs.push_back(i*nocc+i, E[0][i]);
            }
            for (int i = 0;i < nvrt;i++)
            {
                abpairs.push_back(i*nvrt+i, E[0][i+nocc]);
            }

            H.getAB().setDataBySpinAndIrrep({1,1}, {0,0}, abpairs);
            H.getAB().setDataBySpinAndIrrep({0,0}, {0,0}, abpairs);
            H.getIJ().setDataBySpinAndIrrep({1,1}, {0,0}, ijpairs);
            H.getIJ().setDataBySpinAndIrrep({0,0}, {0,0}, ijpairs);

            writeIntegrals( true,  true,  true,  true, H.getABCD());
            writeIntegrals( true,  true,  true, false, H.getABCI());
            writeIntegrals( true,  true, false, false, H.getABIJ());
            writeIntegrals( true, false,  true,  true, H.getAIBC());
            writeIntegrals( true, false,  true, false, H.getAIBJ());
            writeIntegrals( true, false, false, false, H.getAIJK());
            writeIntegrals(false, false,  true,  true, H.getIJAB());
            writeIntegrals(false, false,  true, false, H.getIJAK());
            writeIntegrals(false, false, false, false, H.getIJKL());

            return true;
        }
};

}
}

static const char* spec = R"!(

radius double,
num_electrons int,
num_orbitals int,
dimension? int 3

)!";

REGISTER_TASK(aquarius::jellium::Jellium,"jellium",spec);
