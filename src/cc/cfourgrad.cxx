#include "cfourgrad.hpp"

#include "input/molecule.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

using namespace aquarius::op;
using namespace aquarius::input;
using namespace aquarius::tensor;
using namespace aquarius::task;
using namespace aquarius::time;
using namespace aquarius::symmetry;
using namespace aquarius::integrals;

namespace aquarius
{
namespace cc
{

CFOURGradient::CFOURGradient(const string& name, Config& config)
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.emplace_back("molecule", "molecule");
    reqs.emplace_back(config.get<string>("source")+".D", "D");
    this->addProduct("gradient", "gradient", reqs);
}

bool CFOURGradient::run(task::TaskDAG& dag, const Arena& arena)
{
    auto& molecule = this->template get<Molecule>("molecule");
    auto& D = this->template get<TwoElectronOperator<double>>("D");

    auto& occ = D.occ;
    auto& vrt = D.vrt;
    auto& group = occ.group;

    mkdir(".cfour", 0777);
    chdir(".cfour");

    ofstream ofs("ZMAT");

    ofs << "AQ" << endl;
    for (auto& atom : molecule.getAtoms())
    {
        auto& center = atom.getCenter();
        string sym = toupper(center.getElement().getSymbol());

        for (auto& pos : center.getCenters())
        {
            ofs << printos("%s % 25.18f % 25.18f % 25.18f",
                           sym, pos[0], pos[1], pos[2]) << endl;
        }
    }
    ofs << endl;
    ofs << "*CFOUR(BASIS=SPECIAL" << endl;
    ofs << "CALC=CCSD" << endl;
    ofs << "CC_PROG=NCC" << endl;
    ofs << "COORDS=CARTESIAN" << endl;
    ofs << "UNITS=BOHR)" << endl;
    ofs << endl;
    for (auto& atom : molecule.getAtoms())
    {
        auto& center = atom.getCenter();
        string sym = toupper(center.getElement().getSymbol());

        for (auto& pos : center.getCenters())
        {
            ofs << printos("%s:%p", sym, &atom) << endl;
        }
    }
    ofs << endl;

    ofs.close();

    ofs.open("GENBAS");

    for (auto& atom : molecule.getAtoms())
    {
        auto& center = atom.getCenter();
        string sym = toupper(center.getElement().getSymbol());
        string name = str("%s:%p", sym, &atom);
        vector<Shell> shells(atom.getShellsBegin(), atom.getShellsEnd());

        ofs << name << endl;
        ofs << "AQ" << endl;
        ofs << endl;
        ofs << printos("%3d", shells.size()) << endl;
        for (auto& shell : shells) ofs << printos("%5d", shell.getL());
        ofs << endl;
        for (auto& shell : shells) ofs << printos("%5d", shell.getNContr());
        ofs << endl;
        for (auto& shell : shells) ofs << printos("%5d", shell.getNPrim());
        ofs << endl;
        ofs << endl;
        for (auto& shell : shells)
        {
            auto& exp = shell.getExponents();
            auto& coef = shell.getCoefficients();
            int n = exp.size();
            int m = coef.size()/n;

            for (int i = 0;i < n;)
            {
                for (int j = 0;j < 5 && i < n;i++, j++)
                {
                    ofs << printos("%14.6f", exp[i]);
                }
                ofs << endl;
            }
            ofs << endl;

            for (int k = 0;k < n;k++)
            {
                for (int i = 0;i < m;)
                {
                    for (int j = 0;j < 5 && i < m;i++, j++)
                    {
                        ofs << printos("% 11.7f", coef[k+i*n]);
                    }
                    ofs << endl;
                }
            }
            ofs << endl;
        }
    }

    ofs.close();

    /*
    DIR* dir = opendir(".");
    while (true)
    {
        dirent* ent = readdir(dir);
        if (!ent) break;
        unlink(ent->d_name);
    }
    closedir(dir);
    */
    chdir("..");
    //rmdir(".cfour");

    return true;
}

}
}

static const char* spec = R"!(

source
    string

)!";

REGISTER_TASK(aquarius::cc::CFOURGradient,"cfourgrad",spec);
