#include "cfourgrad.hpp"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>

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

class JOBARC
{
    public:
        JOBARC()
        {
            ja = fopen("JAINDX", "r+b");
            if (!ja)
            {
                perror("fopen: JAINDX");
                abort();
            }

            jo = fopen("JOBARC", "r+b");
            if (!jo)
            {
                perror("fopen: JOBARC");
                abort();
            }

            fseek(ja, 0, SEEK_END);
            size_t len = ftell(ja);
            fseek(ja, 0, SEEK_SET);

            if (len == 2*4 + 8000 + 2001*4)
            {
                recsize = 4;
                intsize = 4;
            }
            else if (len == 2*8 + 8000 + 2001*4)
            {
                recsize = 8;
                intsize = 4;
            }
            else if (len == 2*4 + 8000 + 2001*8)
            {
                recsize = 4;
                intsize = 8;
            }
            else if (len == 2*8 + 8000 + 2001*8)
            {
                recsize = 8;
                intsize = 8;
            }
            else
            {
                abort();
            }

            cout << recsize << " " << intsize << endl;

            size_t reclen = read_int(ja, recsize);
            assert(reclen == 8000 + 2001*intsize);

            for (int i = 0;i < 1000;i++)
            {
                fread(labels[i], 8, 1, ja);
                labels[i][8] = 0;
                cout << labels[i] << endl;
            }

            for (int i = 0;i < 1000;i++)
            {
                off[i] = read_int(ja, intsize);
            }

            for (int i = 0;i < 1000;i++)
            {
                size[i] = read_int(ja, intsize);
            }

            nrecs = read_int(ja, intsize);

            assert(nrecs > 0 && nrecs < 1000);
            assert(strncmp(labels[nrecs  ], "OPENSLOT", 8) == 0);
            assert(strncmp(labels[nrecs-1], "OPENSLOT", 8) != 0);

            reclen = read_int(ja, recsize);
            assert(reclen == 8000 + 2001*intsize);
        }

        ~JOBARC()
        {
            write_int(ja, 8000 + 2001*intsize, recsize);

            for (int i = 0;i < 1000;i++)
            {
                fwrite(labels[i], 8, 1, ja);
            }

            for (int i = 0;i < 1000;i++)
            {
                write_int(ja, off[i], intsize);
            }

            for (int i = 0;i < 1000;i++)
            {
                write_int(ja, size[i], intsize);
            }

            write_int(ja, nrecs, intsize);

            write_int(ja, 8000 + 2001*intsize, recsize);

            fclose(ja);
            fclose(jo);
        }

        void put(const string& label, const string& value)
        {
            size_t size = value.size();
            size_t pad = (size%intsize == 0 ? 0 : intsize-size%intsize);
            string pad_val(size+pad, ' ');
            copy_n(value.begin(), size, pad_val.begin());
            put(label, size+pad, &pad_val[0]);
        }

        template <typename T>
        enable_if_integral_t<T> put(const string& label, const T& value)
        {
            if (intsize == 4)
            {
                int32_t ival = value;
                put(label, 4, &ival);
            }
            else
            {
                int64_t ival = value;
                put(label, 8, &ival);
            }
        }

        template <typename T>
        enable_if_floating_point<T> put(const string& label, const T& value)
        {
            double dval = value;
            put(label, 8, &dval);
        }

        template <typename T>
        enable_if_integral_t<T> put(const string& label, const vector<T>& contents)
        {
            if (intsize == sizeof(T))
            {
                put(label, intsize*contents.size(), contents.data());
            }
            else if (intsize == 4)
            {
                vector<int32_t> icont(contents.begin(), contents.end());
                put(label, 4*icont.size(), icont.data());
            }
            else
            {
                vector<int64_t> icont(contents.begin(), contents.end());
                put(label, 4*icont.size(), icont.data());
            }
        }

        template <typename T>
        enable_if_floating_point_t<T> put(const string& label, const vector<T>& contents)
        {
            if (is_same<T,double>::value)
            {
                put(label, 8*contents.size(), contents.data());
            }
            else
            {
                vector<double> dcont(contents.begin(), contents.end());
                put(label, 8*dcont.size(), dcont.data());
            }
        }

        void put(const string& label, size_t size, const void* contents)
        {
            assert(size > 0 && size%intsize == 0);
            size /= intsize;

            int idx = get_idx(label);

            if (idx == -1)
            {
                idx = get_idx("OPENSLOT");
                assert(idx > 0);

                strncpy(labels[idx], label.c_str(), 8);
                off[idx] = off[idx-1]+this->size[idx-1];
                this->size[idx] = size;
            }

            assert(size <= this->size[idx]);

            fseek(jo, intsize*(off[idx]-1), SEEK_SET);
            fwrite(contents, intsize, size, jo);
        }

        template <typename T>
        T get(const string& label)
        {
            T value;
            get(label, value);
            return value;
        }

        void get(const string& label, string& value)
        {
            int idx = get_idx(label);
            assert(idx != -1);

            value.resize(size[idx]*intsize);
            get(label, size[idx]*intsize, &value[0]);
            while (value.back() == ' ') value.pop_back();
        }

        template <typename T>
        enable_if_integral_t<T> get(const string& label, T& value)
        {
            if (intsize == 4)
            {
                int32_t ival;
                get(label, 4, &ival);
                value = ival;
            }
            else
            {
                int64_t ival;
                get(label, 8, &ival);
                value = ival;
            }
        }

        template <typename T>
        enable_if_floating_point_t<T> get(const string& label, T& value)
        {
            double dval;
            get(label, 8, &dval);
            value = dval;
        }

        template <typename T>
        enable_if_integral_t<T> get(const string& label, vector<T>& contents)
        {
            int idx = get_idx(label);
            assert(idx != -1);

            if (intsize == sizeof(T))
            {
                contents.resize(size[idx]);
                get(label, intsize*contents.size(), contents.data());
            }
            else if (intsize == 4)
            {
                vector<int32_t> icont(size[idx]);
                get(label, 4*icont.size(), icont.data());
                contents.assign(icont.begin(), icont.end());
            }
            else
            {
                vector<int64_t> icont(size[idx]);
                get(label, 8*icont.size(), icont.data());
                contents.assign(icont.begin(), icont.end());
            }
        }

        template <typename T>
        enable_if_floating_point_t<T> get(const string& label, vector<T>& contents)
        {
            int idx = get_idx(label);
            assert(idx != -1);

            if (is_same<double,T>::value)
            {
                contents.resize(size[idx]);
                get(label, 8*contents.size(), contents.data());
            }
            else
            {
                vector<double> dcont(size[idx]);
                get(label, 8*dcont.size(), dcont.data());
                contents.assign(dcont.begin(), dcont.end());
            }
        }

        void get(const string& label, size_t size, void* contents)
        {
            assert(size%intsize == 0);
            size /= intsize;

            int idx = get_idx(label);
            assert(idx != -1);
            assert(size <= this->size[idx]);

            fseek(jo, intsize*(off[idx]-1), SEEK_SET);
            fread(contents, intsize, size, jo);
        }

    protected:
        FILE* ja;
        FILE* jo;
        char labels[1000][9];
        size_t off[1000];
        size_t size[1000];
        int nrecs;
        int recsize;
        int intsize;

        int get_idx(string label)
        {
            while (label.size() < 8) label.push_back(' ');

            for (int idx = 0;idx < 1000;idx++)
            {
                if (strncmp(label.c_str(), labels[idx], 8) == 0) return idx;
            }

            return -1;
        }

        size_t read_int(FILE* fd, int size)
        {
            if (size == 4)
            {
                int32_t x;
                fread(&x, 4, 1, fd);
                return x;
            }
            else
            {
                int64_t x;
                fread(&x, 4, 1, fd);
                return x;
            }
        }

        void write_int(FILE* fd, size_t val, int size)
        {
            if (size == 4)
            {
                int32_t x = val;
                fwrite(&x, 4, 1, fd);
            }
            else
            {
                int64_t x = val;
                fwrite(&x, 4, 1, fd);
            }
        }
};

CFOURGradient::CFOURGradient(const string& name, Config& config)
: Task(name, config), config(config.clone())
{
    vector<Requirement> reqs;
    reqs.emplace_back("molecule", "molecule");
    addProduct("gradient", "gradient", reqs);
    addProduct("moints", "H", reqs);
}

bool CFOURGradient::run(task::TaskDAG& dag, const Arena& arena)
{
    if (getProduct("gradient").getRequirements().size() == 1)
    {
        getProduct("gradient").addRequirement(config.get<string>("source")+".D", "D");
        config.remove("source");

        clean();
        if (mkdir(".cfour", 0777))
        {
            perror("mkdir");
            abort();
        }
        chdir(".cfour");

        writeZMAT();
        writeGENBAS();

        execute(arena, "xjoda");
        execute(arena, "xvmol");
        execute(arena, "xvmol2ja");
        execute(arena, "xvscf");
        execute(arena, "xvtran");
        execute(arena, "xintprc");
        execute(arena, "xint");

        readIntegrals(arena);

        for (auto& c : config.find("*"))
        {
            Task& t = dag.addTask(arena, c.first, getName(), c.second);
        }

        return false;
    }

    writeDensity();

    execute(arena, "xdens");
    execute(arena, "xanti");
    execute(arena, "xbcktrn");
    execute(arena, "xvdint");

    /*
     * Read gradient
     */

    chdir("..");
    //clean();

    return true;
}

void CFOURGradient::writeZMAT()
{
    auto& molecule = get<Molecule>("molecule");

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
    ofs << "CC_PROG=MRCC" << endl;
    ofs << "SCF_CONV=" << lround(-log10(config.get<double>("scf.convergence"))) << endl;
    ofs << "SCF_DAMPING=" << lround(config.get<double>("scf.diis.damping")*1000) << endl;
    ofs << "SCF_EXPSTART=" << config.get<int>("scf.diis.start") << endl;
    ofs << "SCF_EXPORDER=" << config.get<int>("scf.diis.order") << endl;
    ofs << "SCF_MAXCYC=" << config.get<int>("scf.max_iterations") << endl;
    ofs << "FROZEN_CORE=" << (config.get<bool>("scf.frozen_core") ? "ON" : "OFF") << endl;
    ofs << "SYM=OFF" << endl;
    ofs << "REF=UHF" << endl;
    ofs << "DERIV_LEV=1" << endl;
    ofs << "MULT=" << (molecule.getNumAlphaElectrons()-
                       molecule.getNumBetaElectrons()+1) << endl;
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
    config.remove("scf");
}

void CFOURGradient::writeGENBAS()
{
    // TWOPI_N34 = (2 pi)^(-3/4)
    const double TWOPI_N34 = 0.25197943553838073034791409490358;

    auto& molecule = get<Molecule>("molecule");

    ofstream ofs("GENBAS");

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
            int L = shell.getL();

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
                double fac = 1.0/(TWOPI_N34*pow(2*exp[k],0.5*L+0.75));
                for (int i = 0;i < m;)
                {
                    for (int j = 0;j < 5 && i < m;i++, j++)
                    {
                        ofs << printos("% 11.7f", coef[k+i*n]*fac);
                    }
                    ofs << endl;
                }
            }
            ofs << endl;
        }
    }

    ofs.close();
}

void CFOURGradient::execute(const Arena& arena, const string& module)
{
    log(arena) << "executing external program: " << module << endl;

    string command = module + " 2>&1";
    FILE* child = popen(command.c_str(), "r");
    if (!child)
    {
        perror("popen");
        abort();
    }

    if (arena.rank == 0)
    {
        char c;
        while ((c = getc(child)) != EOF) putchar(c);
    }

    pclose(child);
}

void CFOURGradient::readIntegrals(const Arena& arena)
{
    auto& molecule = get<Molecule>("molecule");

    ifstream ifs("fort.55");

    int N, nelec;
    ifs >> N >> nelec;
    int ndrop = (molecule.getNumAlphaElectrons()+molecule.getNumBetaElectrons()-nelec)/2;

    Space occ(PointGroup::C1(), {molecule.getNumAlphaElectrons()-ndrop},
                                {molecule.getNumBetaElectrons()-ndrop});
    Space vrt(PointGroup::C1(), {N-occ.nalpha[0]}, {N-occ.nbeta[0]});

    auto& H = put("H", new TwoElectronOperator<double>("H", arena, occ, vrt));

    for (int i = 0;i < N+1;i++)
    {
        int ignore;
        ifs >> ignore;
    }

    int na = occ.nalpha[0];
    int nb = occ.nbeta[0];

    readIntegrals(ifs, H.getIJKL()({0,2},{0,2})({0,0,0,0}), 0, 0, 0, 0, false, false);

    readIntegrals(ifs, H.getIJAK()({0,2},{1,1})({0,0,0,0}), 0, 0, na, 0, false, true);
    H.getAIJK()({1,1},{0,2})({0,0,0,0})["rspq"] = H.getIJAK()({0,2},{1,1})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABIJ()({2,0},{0,2})({0,0,0,0}), na, na, 0, 0, false, false);
    H.getIJAB()({0,2},{2,0})({0,0,0,0})["rspq"] = H.getABIJ()({2,0},{0,2})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getAIBJ()({1,1},{1,1})({0,0,0,0}), na, 0, na, 0, true, true);

    readIntegrals(ifs, H.getABCI()({2,0},{1,1})({0,0,0,0}), na, na, na, 0, false, false);
    H.getAIBC()({1,1},{2,0})({0,0,0,0})["rspq"] = H.getABCI()({2,0},{1,1})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABCD()({2,0},{2,0})({0,0,0,0}), na, na, na, na, false, false);

    if (arena.rank == 0)
    {
        double value;
        int p, q, r, s;
        ifs >> value >> p >> q >> r >> s;
        assert(p == 0 && q == 0 && r == 0 && s == 0);
    }

    readIntegrals(ifs, H.getIJKL()({0,0},{0,0})({0,0,0,0}), 0, 0, 0, 0, false, false);

    readIntegrals(ifs, H.getIJAK()({0,0},{0,0})({0,0,0,0}), 0, 0, nb, 0, false, true);
    H.getAIJK()({0,0},{0,0})({0,0,0,0})["rspq"] = H.getIJAK()({0,0},{0,0})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABIJ()({0,0},{0,0})({0,0,0,0}), nb, nb, 0, 0, false, false);
    H.getIJAB()({0,0},{0,0})({0,0,0,0})["rspq"] = H.getABIJ()({0,0},{0,0})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getAIBJ()({0,0},{0,0})({0,0,0,0}), nb, 0, nb, 0, true, true);

    readIntegrals(ifs, H.getABCI()({0,0},{0,0})({0,0,0,0}), nb, nb, nb, 0, false, false);
    H.getAIBC()({0,0},{0,0})({0,0,0,0})["rspq"] = H.getABCI()({0,0},{0,0})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABCD()({0,0},{0,0})({0,0,0,0}), nb, nb, nb, nb, false, false);

    if (arena.rank == 0)
    {
        double value;
        int p, q, r, s;
        ifs >> value >> p >> q >> r >> s;
        assert(p == 0 && q == 0 && r == 0 && s == 0);
    }

    readIntegrals(ifs, H.getIJKL()({0,1},{0,1})({0,0,0,0}), 0, 0, 0, 0, false, false);

    readIntegrals(ifs, H.getIJAK()({0,1},{0,1})({0,0,0,0}), 0, 0, nb, 0, false, true);
    H.getAIJK()({0,1},{0,1})({0,0,0,0})["rspq"] = H.getIJAK()({0,1},{0,1})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getIJAK()({0,1},{1,0})({0,0,0,0}), 0, 0, na, 0, false, false);
    H.getAIJK()({1,0},{0,1})({0,0,0,0})["rspq"] = H.getIJAK()({0,1},{1,0})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABIJ()({1,0},{0,1})({0,0,0,0}), na, nb, 0, 0, false, false);
    H.getIJAB()({0,1},{1,0})({0,0,0,0})["rspq"] = H.getABIJ()({1,0},{0,1})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getAIBJ()({1,0},{1,0})({0,0,0,0}), na, 0, na, 0, false, false);
    readIntegrals(ifs, H.getAIBJ()({0,1},{0,1})({0,0,0,0}), nb, 0, nb, 0, true, true);

    readIntegrals(ifs, H.getABCI()({1,0},{1,0})({0,0,0,0}), na, nb, na, 0, false, false);
    H.getAIBC()({1,0},{1,0})({0,0,0,0})["rspq"] = H.getABCI()({1,0},{1,0})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABCI()({1,0},{0,1})({0,0,0,0}), na, nb, nb, 0, false, true);
    H.getAIBC()({0,1},{1,0})({0,0,0,0})["rspq"] = H.getABCI()({1,0},{0,1})({0,0,0,0})["pqrs"];

    readIntegrals(ifs, H.getABCD()({1,0},{1,0})({0,0,0,0}), na, nb, na, nb, false, false);

    H.getAIBJ()({1,0},{0,1})["AibJ"] = -H.getABIJ()({1,0},{0,1})["AbJi"];
    H.getAIBJ()({0,1},{1,0})["aIBj"] = -H.getABIJ()({1,0},{0,1})["BaIj"];

    if (arena.rank == 0)
    {
        double value;
        int p, q, r, s;
        ifs >> value >> p >> q >> r >> s;
        assert(p == 0 && q == 0 && r == 0 && s == 0);
    }

    readIntegrals(ifs, H.getIJ()({0,1},{0,1})({0,0}), 0, 0);
    readIntegrals(ifs, H.getAB()({1,0},{1,0})({0,0}), na, na);
    readIntegrals(ifs, H.getIA()({0,1},{1,0})({0,0}), 0, na);
    H.getAI()({1,0},{0,1})({0,0})["qp"] = H.getIA()({0,1},{1,0})({0,0})["pq"];

    if (arena.rank == 0)
    {
        double value;
        int p, q, r, s;
        ifs >> value >> p >> q >> r >> s;
        assert(p == 0 && q == 0 && r == 0 && s == 0);
    }

    readIntegrals(ifs, H.getIJ()({0,0},{0,0})({0,0}), 0, 0);
    readIntegrals(ifs, H.getAB()({0,0},{0,0})({0,0}), nb, nb);
    readIntegrals(ifs, H.getIA()({0,0},{0,0})({0,0}), 0, nb);
    H.getAI()({0,0},{0,0})({0,0})["qp"] = H.getIA()({0,0},{0,0})({0,0})["pq"];

    SpinorbitalTensor<double> D("D", arena, PointGroup::C1(), {vrt,occ}, {0,1}, {0,1});

    if (arena.rank == 0)
    {
        vector<kv_pair> pairsa, pairsb;
        for (int i = 0;i < na;i++) pairsa.emplace_back(i+i*na, 1.0);
        for (int i = 0;i < nb;i++) pairsb.emplace_back(i+i*nb, 1.0);
        D({0,1},{0,1})({0,0}).writeRemoteData(pairsa);
        D({0,0},{0,0})({0,0}).writeRemoteData(pairsb);
    }
    else
    {
        D({0,1},{0,1})({0,0}).writeRemoteData();
        D({0,0},{0,0})({0,0}).writeRemoteData();
    }

    H.getAB()["ab"] += H.getAIBJ()["ambn"]*D["nm"];
    H.getAI()["ai"] += H.getAIJK()["amin"]*D["nm"];
    H.getIA()["ia"] += H.getIJAK()["iman"]*D["nm"];
    H.getIJ()["ij"] += H.getIJKL()["imjn"]*D["nm"];

    H.getAI() = 0;
    H.getIA() = 0;

    //this->log(arena) << "ABCD: " << setprecision(15) << H.getABCD()({2,0},{2,0}).norm(2) << endl;
    //this->log(arena) << "AbCd: " << setprecision(15) << H.getABCD()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "abcd: " << setprecision(15) << H.getABCD()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "ABCI: " << setprecision(15) << H.getABCI()({2,0},{1,1}).norm(2) << endl;
    //this->log(arena) << "AbCi: " << setprecision(15) << H.getABCI()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "AbcI: " << setprecision(15) << H.getABCI()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "abci: " << setprecision(15) << H.getABCI()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIBC: " << setprecision(15) << H.getAIBC()({1,1},{2,0}).norm(2) << endl;
    //this->log(arena) << "AiBc: " << setprecision(15) << H.getAIBC()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "aIBc: " << setprecision(15) << H.getAIBC()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "aibc: " << setprecision(15) << H.getAIBC()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "ABIJ: " << setprecision(15) << H.getABIJ()({2,0},{0,2}).norm(2) << endl;
    //this->log(arena) << "AbIj: " << setprecision(15) << H.getABIJ()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "abij: " << setprecision(15) << H.getABIJ()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIBJ: " << setprecision(15) << H.getAIBJ()({1,1},{1,1}).norm(2) << endl;
    //this->log(arena) << "AiBj: " << setprecision(15) << H.getAIBJ()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "aIbJ: " << setprecision(15) << H.getAIBJ()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "AibJ: " << setprecision(15) << H.getAIBJ()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "aIBj: " << setprecision(15) << H.getAIBJ()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "aibj: " << setprecision(15) << H.getAIBJ()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJAB: " << setprecision(15) << H.getIJAB()({0,2},{2,0}).norm(2) << endl;
    //this->log(arena) << "IjAb: " << setprecision(15) << H.getIJAB()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "ijab: " << setprecision(15) << H.getIJAB()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AIJK: " << setprecision(15) << H.getAIJK()({1,1},{0,2}).norm(2) << endl;
    //this->log(arena) << "AiJk: " << setprecision(15) << H.getAIJK()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "aIJk: " << setprecision(15) << H.getAIJK()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "aijk: " << setprecision(15) << H.getAIJK()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJAK: " << setprecision(15) << H.getIJAK()({0,2},{1,1}).norm(2) << endl;
    //this->log(arena) << "IjAk: " << setprecision(15) << H.getIJAK()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "IjaK: " << setprecision(15) << H.getIJAK()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ijak: " << setprecision(15) << H.getIJAK()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJKL: " << setprecision(15) << H.getIJKL()({0,2},{0,2}).norm(2) << endl;
    //this->log(arena) << "IjKl: " << setprecision(15) << H.getIJKL()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ijkl: " << setprecision(15) << H.getIJKL()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AB:   " << setprecision(15) << H.getAB()({1,0},{1,0}).norm(2) << endl;
    //this->log(arena) << "ab:   " << setprecision(15) << H.getAB()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "AI:   " << setprecision(15) << H.getAI()({1,0},{0,1}).norm(2) << endl;
    //this->log(arena) << "ai:   " << setprecision(15) << H.getAI()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IA:   " << setprecision(15) << H.getIA()({0,1},{1,0}).norm(2) << endl;
    //this->log(arena) << "ia:   " << setprecision(15) << H.getIA()({0,0},{0,0}).norm(2) << endl;
    //this->log(arena) << "IJ:   " << setprecision(15) << H.getIJ()({0,1},{0,1}).norm(2) << endl;
    //this->log(arena) << "ij:   " << setprecision(15) << H.getIJ()({0,0},{0,0}).norm(2) << endl;
}

void CFOURGradient::readIntegrals(ifstream& ifs, CTFTensor<double>& H,
                                  int offp, int offq, int offr, int offs,
                                  bool transpq, bool transrs)
{
    auto& len = H.getLengths();
    bool sympq = H.getSymmetry()[0] == AS;
    bool symrs = H.getSymmetry()[2] == AS;
    int64_t numpq = (sympq ? int64_t(len[0])*(len[1]-1)/2 : int64_t(len[0])*len[1]);
    int64_t numrs = (symrs ? int64_t(len[2])*(len[3]-1)/2 : int64_t(len[2])*len[3]);

    if (H.arena.rank == 0)
    {
        vector<kv_pair> pairs(numpq*numrs);

        for (auto& pair : pairs)
        {
            int p, q, r, s;
            double value;
            ifs >> value >> p >> r >> q >> s;
            pair.d = value;

            if (transpq)
            {
                pair.d = -pair.d;
                if (!sympq) swap(p, q);
            }

            if (transrs)
            {
                pair.d = -pair.d;
                if (!symrs) swap(r, s);
            }

            p = p-offp-1;
            q = q-offq-1;
            r = r-offr-1;
            s = s-offs-1;

            assert(p >= 0 && p < len[0]);
            assert(q >= 0 && q < len[1]);
            assert(r >= 0 && r < len[2]);
            assert(s >= 0 && s < len[3]);
            if (sympq) assert(p < q);
            if (symrs) assert(r < s);
            pair.k = ((int64_t(s)*len[2]+r)*len[1]+q)*len[0]+p;
        }

        H.writeRemoteData(pairs);
    }
    else
    {
        H.writeRemoteData();
    }
}

void CFOURGradient::readIntegrals(ifstream& ifs, CTFTensor<double>& H,
                                  int offp, int offq)
{
    auto& len = H.getLengths();

    if (H.arena.rank == 0)
    {
        vector<kv_pair> pairs(len[0]*len[1]);

        for (auto& pair : pairs)
        {
            int p, q, r, s;
            double value;
            ifs >> value >> p >> q >> r >> s;
            p = p-offp-1;
            q = q-offq-1;
            pair.d = value;

            assert(p >= 0 && p < len[0]);
            assert(q >= 0 && q < len[1]);
            assert(r == 0);
            assert(s == 0);
            pair.k = int64_t(q)*len[0]+p;
        }

        H.writeRemoteData(pairs);
    }
    else
    {
        H.writeRemoteData();
    }
}

void CFOURGradient::writeDensity()
{
    auto& D = get<TwoElectronOperator<double>>("D");
    auto& arena = D.arena;

    ofstream ofs;

    int na = D.occ.nalpha[0];
    int nb = D.occ.nbeta[0];

    if (arena.rank == 0) ofs.open("CCDENSITIES");

    if (arena.rank == 0) ofs << " G(IJ,KL)" << endl;
    writeDensity(ofs, D.getIJKL()({0,2},{0,2})({0,0,0,0}), 0, 0, 0, 0, false, false);

    if (arena.rank == 0) ofs << " G(ij,kl)" << endl;
    writeDensity(ofs, D.getIJKL()({0,0},{0,0})({0,0,0,0}), 0, 0, 0, 0, false, false);

    if (arena.rank == 0) ofs << " G(Ij,Kl)" << endl;
    writeDensity(ofs, D.getIJKL()({0,1},{0,1})({0,0,0,0}), 0, 0, 0, 0, false, false);

    {
        CTFTensor<double> G(D.getIJAK()({0,2},{1,1})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIJK()({1,1},{0,2})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(IJ,KA)" << endl;
        writeDensity(ofs, G, 0, 0, na, 0, true, true);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,0},{0,0})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIJK()({0,0},{0,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(ij,ka)" << endl;
        writeDensity(ofs, G, 0, 0, nb, 0, true, true);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,1},{1,0})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIJK()({1,0},{0,1})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(Ij,Ak)" << endl;
        writeDensity(ofs, G, 0, 0, na, 0, false, false);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,1},{0,1})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIJK()({0,1},{0,1})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(Ij,Ka)" << endl;
        writeDensity(ofs, G, 0, 0, nb, 0, false, true);
    }

    {
        CTFTensor<double> G(D.getABIJ()({2,0},{0,2})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getIJAB()({0,2},{2,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(AB,IJ)" << endl;
        writeDensity(ofs, G, na, na, 0, 0, false, false);
    }

    {
        CTFTensor<double> G(D.getABIJ()({0,0},{0,0})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getIJAB()({0,0},{0,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(ab,ij)" << endl;
        writeDensity(ofs, G, nb, nb, 0, 0, false, false);
    }

    {
        CTFTensor<double> G(D.getABIJ()({1,0},{0,1})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getIJAB()({0,1},{1,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(Ab,Ij)" << endl;
        writeDensity(ofs, G, na, nb, 0, 0, false, false);
    }

    if (arena.rank == 0) ofs << " G(AI,BJ)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,1},{1,1})({0,0,0,0}), na, 0, na, 0, false, false);

    if (arena.rank == 0) ofs << " G(ai,bj)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,0},{0,0})({0,0,0,0}), nb, 0, nb, 0, false, false);

    if (arena.rank == 0) ofs << " G(Ai,Bj)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,0},{1,0})({0,0,0,0}), na, 0, na, 0, false, false);

    if (arena.rank == 0) ofs << " G(aI,bJ)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,1},{0,1})({0,0,0,0}), nb, 0, nb, 0, false, false);

    if (arena.rank == 0) ofs << " G(Aj,Ib)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,0},{0,1})({0,0,0,0}), na, 0, nb, 0, false, true);

    if (arena.rank == 0) ofs << " G(aJ,iB)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,1},{1,0})({0,0,0,0}), nb, 0, na, 0, false, true);

    {
        CTFTensor<double> G(D.getABCI()({2,0},{1,1})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIBC()({1,1},{2,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(AB,CI)" << endl;
        writeDensity(ofs, G, na, na, na, 0, false, false);
    }

    {
        CTFTensor<double> G(D.getABCI()({0,0},{0,0})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIBC()({0,0},{0,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(ab,ci)" << endl;
        writeDensity(ofs, G, nb, nb, nb, 0, false, false);
    }

    {
        CTFTensor<double> G(D.getABCI()({1,0},{0,1})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIBC()({0,1},{1,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(Ab,Ic)" << endl;
        writeDensity(ofs, G, na, nb, nb, 0, false, true);
    }

    {
        CTFTensor<double> G(D.getABCI()({1,0},{1,0})({0,0,0,0}));
        0.5*G["pqrs"] += 0.5*D.getAIBC()({1,0},{1,0})({0,0,0,0})["rspq"];
        if (arena.rank == 0) ofs << " G(Ab,Ci)" << endl;
        writeDensity(ofs, G, na, nb, na, 0, false, false);
    }

    if (arena.rank == 0) ofs << " G(AB,CD)" << endl;
    writeDensity(ofs, D.getABCD()({2,0},{2,0})({0,0,0,0}), na, na, na, na, false, false);

    if (arena.rank == 0) ofs << " G(ab,cd)" << endl;
    writeDensity(ofs, D.getABCD()({0,0},{0,0})({0,0,0,0}), nb, nb, nb, nb, false, false);

    if (arena.rank == 0) ofs << " G(Ab,Cd)" << endl;
    writeDensity(ofs, D.getABCD()({1,0},{1,0})({0,0,0,0}), na, nb, na, nb, false, false);

    if (arena.rank == 0) ofs << " D(I,J)  " << endl;
    writeDensity(ofs, D.getIJ()({0,1},{0,1})({0,0}), 0, 0);

    if (arena.rank == 0) ofs << " D(A,B)  " << endl;
    writeDensity(ofs, D.getAB()({1,0},{1,0})({0,0}), na, na);

    {
        CTFTensor<double> G(D.getAI()({1,0},{0,1})({0,0}));
        0.5*G["pq"] += 0.5*D.getIA()({0,1},{1,0})({0,0})["qp"];
        if (arena.rank == 0) ofs << " D(A,I)  " << endl;
        writeDensity(ofs, G, na, 0);
    }

    if (arena.rank == 0) ofs << " D(i,j)  " << endl;
    writeDensity(ofs, D.getIJ()({0,0},{0,0})({0,0}), 0, 0);

    if (arena.rank == 0) ofs << " D(a,b)  " << endl;
    writeDensity(ofs, D.getAB()({0,0},{0,0})({0,0}), nb, nb);

    {
        CTFTensor<double> G(D.getAI()({0,0},{0,0})({0,0}));
        0.5*G["pq"] += 0.5*D.getIA()({0,0},{0,0})({0,0})["qp"];
        if (arena.rank == 0) ofs << " D(a,i)  " << endl;
        writeDensity(ofs, G, nb, 0);
    }
}

void CFOURGradient::writeDensity(ofstream& ofs, const CTFTensor<double>& D,
                                 int offp, int offq, int offr, int offs,
                                 bool transpq, bool transrs)
{
    auto& len = D.getLengths();
    bool sympq = D.getSymmetry()[0] == AS;
    bool symrs = D.getSymmetry()[2] == AS;

    if (D.arena.rank == 0)
    {
        vector<double> pairs;
        D.getAllData(pairs, 0);

        ofs << printos("%20d", pairs.size()) << endl;

        int64_t off = 0;
        for (int ss = 0;ss < len[3];ss++)
        {
            int maxr = (symrs ? ss : len[2]);
            for (int rr = 0;rr < maxr;rr++)
            {
                for (int qq = 0;qq < len[1];qq++)
                {
                    int maxp = (sympq ? qq : len[0]);
                    for (int pp = 0;pp < maxp;pp++)
                    {
                        int p = pp+offp+1;
                        int q = qq+offq+1;
                        int r = rr+offr+1;
                        int s = ss+offs+1;
                        double value = pairs[off++];

                        if (transpq)
                        {
                            value = -value;
                            if (!sympq) swap(p, q);
                        }

                        if (transrs)
                        {
                            value = -value;
                            if (!symrs) swap(r, s);
                        }

                        ofs << printos("%28.20e%4d%4d%4d%4d", value, p, q, r, s) << endl;
                    }
                }
            }
        }
    }
    else
    {
        D.getAllData(0);
    }
}

void CFOURGradient::writeDensity(ofstream& ofs, const CTFTensor<double>& D,
                                 int offp, int offq)
{
    auto& len = D.getLengths();

    if (D.arena.rank == 0)
    {
        vector<double> pairs;
        D.getAllData(pairs, 0);

        ofs << printos("%20d", pairs.size()) << endl;

        int64_t off = 0;
        for (int q = 0;q < len[1];q++)
        {
            for (int p = 0;p < len[0];p++)
            {
                ofs << printos("%28.20e%4d%4d", pairs[off++], p+offp+1, q+offq+1) << endl;
            }
        }
    }
    else
    {
        D.getAllData(0);
    }
}

void CFOURGradient::clean()
{
    system("rm -rf .cfour");
}

}
}

static const char* spec = R"!(

source
    string,
scf?
{
    frozen_core?
        bool false,
    convergence?
        double 1e-7,
    max_iterations?
        int 150,
    conv_type?
        enum { MAXE, RMSE, MAE },
    diis?
    {
        damping?
            double 0.0,
        start?
            int 8,
        order?
            int 6,
        jacobi?
            bool false
    }
},
*+

)!";

REGISTER_TASK(aquarius::cc::CFOURGradient,"cfourgrad",spec);
