#include "../../frameworks/input/molecule.hpp"
#include "../../frameworks/molecule/molecule.hpp"

#include "../../frameworks/input/config.hpp"
#include "../../frameworks/util/global.hpp"
using namespace aquarius::input;

#include "../../frameworks/task/task.hpp"
using namespace aquarius::task;

namespace aquarius
{

class MoleculeTask : public task::Task
{
    protected:
        Config config;

    public:
        MoleculeTask(const string& name, input::Config& config)
        : Task(name, config), config(config)
        {
            addProduct(Product("molecule", "molecule"));
        }

        bool run(task::TaskDAG& dag, const Arena& arena)
        {
            put("molecule", new Molecule(config, arena));
            return true;
        }
};

}

static const char* spec = R"(

angstrom2bohr?
    double 1.88972612456506198632428439,
units?
    enum { angstrom, bohr },
coords?
    enum { internal, cartesian },
multiplicity?
    int 1,
charge?
    int 0,
subgroup?
    enum
    {
        C1, full,
        Cs, Ci, C2, C2v, C2h, D2, D2h,
        C3, C4, C5, C6, C3v, C4v, C5v, C6v,
        C3h, C4h, C5h, C6h, D3, D4, D5, D6,
        D3h, D4h, D5h, D6h, S4, S6, Td, Oh, Ih
    },
atom+
{
    basis_set? string,
    truncation? string,
    # ZMAT specification, e.g. H 2 R 1 A
    # or xyz position
    *+
},
basis?
{
    contaminants?
        bool false,
    spherical?
        bool true,
    basis_set? string,
    truncation*
    {
        # elements affected, e.g. H-Ne
        *,
        # same format as molecule.atom.truncation
        *
    }
}

)";

REGISTER_TASK(aquarius::MoleculeTask,"molecule",spec);
