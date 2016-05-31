#include "cfourgrad.hpp"

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

            fclose(fd);
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
: Task(name, config)
{
    vector<Requirement> reqs;
    reqs.emplace_back("molecule", "molecule");
    addProduct("gradient", "gradient", reqs);
}

bool CFOURGradient::run(task::TaskDAG& dag, const Arena& arena)
{
    /*
     * Create additional dependencies (density) or skip to below if already
     * satisfied.
     */

    if (mkdir(".cfour", 0777))
    {
        perror("mkdir");
        abort();
    }
    chdir(".cfour");

    writeZMAT(molecule);
    writeGENBAS(molecule);

    execute(arena, "xjoda");
    execute(arena, "xvmol");
    execute(arena, "xvmol2ja");
    execute(arena, "xvscf");
    execute(arena, "xvtran");
    execute(arena, "xintprc");

    /*
     * Read integrals
     */

    /*
     * Create CC task and connect dependencies
     */

    /*
     * Yield and run CC
     */

    writeDensity(D);

    execute(arena, "xdens");
    execute(arena, "xanti");
    execute(arena, "xbcktrn");
    execute(arena, "xvdint");

    /*
     * Read gradient
     */

    clean();

    chdir("..");
    if (rmdir(".cfour"))
    {
        perror("rmdir");
        abort();
    }

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
    ofs << "SYM=OFF" << endl;
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

    Space occ(PointGroup::C1, {molecule.getNumAlphaElectrons()-ndrop},
                              {molecule.getNumBetaElectrons()-ndrop});
    Space vrt(PointGroup::C1, {N-occ.nalpha[0]}, {N-occ.nbeta[0]});

    auto& H = put("H", new TwoElectronOperator<double>("H", arena, occ, vrt));

    for (int i = 0;i < N+1;i++)
    {
        int ignore;
        ifs >> ignore;
    }


}

void CFOURGradient::readIntegrals(ifstream& ifs, CTFTensor<double>& H, bool transpq, bool transrs)
{
    auto& len = H.getLengths();
    bool sympq = H.getSymmetry()[0] == AS;
    bool symrs = H.getSymmetry()[2] == AS;

    if (H.arena.rank == 0)
    {
        vector<kv_pair> pairs;

        ofs << printos("%20d", pairs.size()) << endl;

        for (auto& pair : pairs)
        {
            int p = pair.k%len[0];
            pair.k /= len[0];
            int q = pair.k%len[1];
            pair.k /= len[1];
            int r = pair.k%len[2];
            pair.k /= len[2];
            int s = pair.k;

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

            ofs << printos("%28.20e%4d%4d%4d%4d", pair.d, p, q, r, s) << endl;
        }
    }
    else
    {
        D.getAllData(0);
    }
}

void CFOURGradient::writeDensity()
{
    auto& D = get<TwoElectronOperator<double>>("D");
    auto& arena = D.arena;

    ofstream ofs("CCDENSITIES");

    if (arena.rank == 0) ofs << " G(IJ,KL)" << endl;
    writeDensity(ofs, D.getIJKL()({0,2},{0,2})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(ij,kl)" << endl;
    writeDensity(ofs, D.getIJKL()({0,0},{0,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(Ij,Kl)" << endl;
    writeDensity(ofs, D.getIJKL()({0,1},{0,1})({0,0,0,0}), false, false);

    {
        CTFTensor<double> G(D.getIJAK()({0,2},{1,1})({0,0,0,0}));
        G["pqrs"] += D.getAIJK()({1,1},{0,2})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(IJ,KA)" << endl;
        writeDensity(ofs, G, true, true);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,0},{0,0})({0,0,0,0}));
        G["pqrs"] += D.getAIJK()({0,0},{0,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(ij,ka)" << endl;
        writeDensity(ofs, G, true, true);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,1},{1,0})({0,0,0,0}));
        G["pqrs"] += D.getAIJK()({1,0},{0,1})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(Ij,Ak)" << endl;
        writeDensity(ofs, G, false, false);
    }

    {
        CTFTensor<double> G(D.getIJAK()({0,1},{0,1})({0,0,0,0}));
        G["pqrs"] += D.getAIJK()({0,1},{0,1})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(Ij,Ka)" << endl;
        writeDensity(ofs, G, false, true);
    }

    {
        CTFTensor<double> G(D.getABIJ()({2,0},{0,2})({0,0,0,0}));
        G["pqrs"] += D.getIJAB()({0,2},{2,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(AB,IJ)" << endl;
        writeDensity(ofs, G, false, false);
    }

    {
        CTFTensor<double> G(D.getABIJ()({0,0},{0,0})({0,0,0,0}));
        G["pqrs"] += D.getIJAB()({0,0},{0,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(ab,ij)" << endl;
        writeDensity(ofs, G, false, false);
    }

    {
        CTFTensor<double> G(D.getABIJ()({1,0},{0,1})({0,0,0,0}));
        G["pqrs"] += D.getIJAB()({0,1},{1,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(Ab,Ij)" << endl;
        writeDensity(ofs, G, false, false);
    }

    if (arena.rank == 0) ofs << " G(AI,BJ)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,1},{1,1})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(ai,bj)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,0},{0,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(Ai,Bj)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,0},{1,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(aI,bJ)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,1},{0,1})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(Aj,Ib)" << endl;
    writeDensity(ofs, D.getAIBJ()({1,0},{0,1})({0,0,0,0}), false, true);

    if (arena.rank == 0) ofs << " G(aJ,iB)" << endl;
    writeDensity(ofs, D.getAIBJ()({0,1},{1,0})({0,0,0,0}), false, true);

    {
        CTFTensor<double> G(D.getABCI()({2,0},{1,1})({0,0,0,0}));
        G["pqrs"] += D.getAIBC()({1,1},{2,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(AB,CI)" << endl;
        writeDensity(ofs, G, false, false);
    }

    {
        CTFTensor<double> G(D.getABCI()({0,0},{0,0})({0,0,0,0}));
        G["pqrs"] += D.getAIBC()({0,0},{0,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(ab,ci)" << endl;
        writeDensity(ofs, G, false, false);
    }

    {
        CTFTensor<double> G(D.getABCI()({1,0},{0,1})({0,0,0,0}));
        G["pqrs"] += D.getAIBC()({0,1},{1,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(Ab,Ic)" << endl;
        writeDensity(ofs, G, false, true);
    }

    {
        CTFTensor<double> G(D.getABCI()({1,0},{1,0})({0,0,0,0}));
        G["pqrs"] += D.getAIBC()({1,0},{1,0})({0,0,0,0})["srpq"];
        if (arena.rank == 0) ofs << " G(Ab,Ci)" << endl;
        writeDensity(ofs, G, false, false);
    }

    if (arena.rank == 0) ofs << " G(AB,CD)" << endl;
    writeDensity(ofs, D.getABCD()({2,0},{2,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(ab,cd)" << endl;
    writeDensity(ofs, D.getABCD()({0,0},{0,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " G(Ab,Cd)" << endl;
    writeDensity(ofs, D.getABCD()({1,0},{1,0})({0,0,0,0}), false, false);

    if (arena.rank == 0) ofs << " D(I,J)  " << endl;
    writeDensity(ofs, D.getIJ()({0,1},{0,1})({0,0}));

    if (arena.rank == 0) ofs << " D(A,B)  " << endl;
    writeDensity(ofs, D.getAB()({1,0},{1,0})({0,0}));

    {
        CTFTensor<double> G(D.getAI()({1,0},{0,1})({0,0}));
        G["pq"] += D.getIA()({0,1},{1,0})({0,0}))["qp"];
        if (arena.rank == 0) ofs << " D(A,I)  " << endl;
        writeDensity(ofs, G);
    }

    if (arena.rank == 0) ofs << " D(i,j)  " << endl;
    writeDensity(ofs, D.getIJ()({0,0},{0,0})({0,0}));

    if (arena.rank == 0) ofs << " D(a,b)  " << endl;
    writeDensity(ofs, D.getAB()({0,0},{0,0})({0,0}));

    {
        CTFTensor<double> G(D.getAI()({0,0},{0,0})({0,0}));
        G["pq"] += D.getIA()({0,0},{0,0})({0,0}))["qp"];
        if (arena.rank == 0) ofs << " D(a,i)  " << endl;
        writeDensity(ofs, G);
    }
}

void CFOURGradient::writeDensity(ofstream& ofs, const CTFTensor<double>& D, bool transpq, bool transrs)
{
    auto& len = D.getLengths();
    bool sympq = D.getSymmetry()[0] == AS;
    bool symrs = D.getSymmetry()[2] == AS;

    if (D.arena.rank == 0)
    {
        vector<kv_pair> pairs;
        D.getAllData(pairs, 0);

        ofs << printos("%20d", pairs.size()) << endl;

        for (auto& pair : pairs)
        {
            int p = pair.k%len[0];
            pair.k /= len[0];
            int q = pair.k%len[1];
            pair.k /= len[1];
            int r = pair.k%len[2];
            pair.k /= len[2];
            int s = pair.k;

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

            ofs << printos("%28.20e%4d%4d%4d%4d", pair.d, p, q, r, s) << endl;
        }
    }
    else
    {
        D.getAllData(0);
    }
}

void CFOURGradient::writeDensity(ofstream& ofs, const CTFTensor<double>& D)
{
    auto& len = D.getLengths();

    if (D.arena.rank == 0)
    {
        vector<kv_pair> pairs;
        D.getAllData(pairs, 0);

        ofs << printos("%20d", pairs.size()) << endl;

        for (auto& pair : pairs)
        {
            int p = pair.k%len[0];
            int q = pair.k/len[0];

            ofs << printos("%28.20e%4d%4d", pair.d, p, q) << endl;
        }
    }
    else
    {
        D.getAllData(0);
    }
}

void CFOURGradient::clean()
{
    DIR* dir = opendir(".");
    if (!dir)
    {
        perror("opendir");
        abort();
    }
    while (true)
    {
        dirent* ent = readdir(dir);
        if (!ent) break;
        struct stat s;
        lstat(ent->d_name, &s);
        if (S_ISREG(s.st_mode))
        {
            if (unlink(ent->d_name))
            {
                perror("unlink");
                abort();
            }
        }
    }
    closedir(dir);
}

}
}

static const char* spec = R"!(

source
    string

)!";

REGISTER_TASK(aquarius::cc::CFOURGradient,"cfourgrad",spec);
