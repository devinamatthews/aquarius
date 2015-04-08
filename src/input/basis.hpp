#ifndef _AQUARIUS_INPUT_BASIS_HPP_
#define _AQUARIUS_INPUT_BASIS_HPP_

#include "util/global.hpp"

#include "integrals/shell.hpp"

#include "molecule.hpp"

namespace aquarius
{
namespace input
{

class BasisSetNotFoundError : public runtime_error
{
    public:
        BasisSetNotFoundError(const string& element) : runtime_error(element) {}
};

class BasisSetFormatError : public runtime_error
{
    private:
        static string buildString(const string& file, const string& what_arg, const int lineno)
        {
            ostringstream os;
            os << file << ": " << what_arg << ": line " << lineno;
            return os.str();
        }

    public:
        BasisSetFormatError(const string& file, const string& what_arg, const int lineno)
        : runtime_error(buildString(file, what_arg, lineno)) {}
};

class BasisSet
{
    private:
        struct ShellBasis
        {
            vector<double> exponents;
            vector<double> coefficients;
            int nprim;
            int ncontr;
            int L;
        };

        map<string,vector<ShellBasis>> atomBases;

        void readBasisSet(const string& file);

        template <typename T>
        T readValue(istream& is, const string& file, int& lineno)
        {
            string line;
            istringstream iss;
            char c;
            T v;

            for (int j = 0;j < 1;)
            {
                line = readLine(is, file, lineno);
                iss.str(line); iss.clear();
                for (;j < 1;j++)
                {
                    do
                    {
                        c = iss.peek();
                        if (c != ' ' && c != '\t') break;
                        iss.get();
                    } while (true);

                    if (iss.eof()) break;

                    if (!(iss >> v)) throw BasisSetFormatError(file, "parse error", lineno);
                }
            }

            return v;
        }

        template <typename T>
        vector<T> readValues(istream& is, const string& file, int& lineno, int n)
        {
            string line;
            istringstream iss;
            char c;
            vector<T> v(n);

            for (int j = 0;j < n;)
            {
                line = readLine(is, file, lineno);
                iss.str(line); iss.clear();
                for (;j < n;j++)
                {
                    do
                    {
                        c = iss.peek();
                        if (c != ' ' && c != '\t') break;
                        iss.get();
                    } while (true);

                    if (iss.eof()) break;

                    if (!(iss >> v[j])) throw BasisSetFormatError(file, "parse error", lineno);
                }
            }

            return v;
        }

        string readLine(istream& is, const string& file, int& lineno);

    public:
        BasisSet() {}

        BasisSet(const string& file);

        void apply(Atom& atom, bool spherical = true, bool contaminants = false);

        void apply(Molecule& molecule, bool spherical = true, bool contaminants = false);
};

}
}

#endif
