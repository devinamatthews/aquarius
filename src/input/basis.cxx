#include "basis.hpp"

using namespace aquarius::integrals;

namespace aquarius
{
namespace input
{

BasisSet::BasisSet(const string& file)
{
    readBasisSet(file);
}

void BasisSet::readBasisSet(const string& file)
{
    ifstream ifs(file.c_str());

    if (!ifs) throw BasisSetNotFoundError(file);

    string line;
    for (int lineno = 1;getline(ifs, line);lineno++)
    {
        // skip blank and comment lines
        if (line.find_first_not_of(" \t\n") == string::npos || line[0] == '!') continue;

        size_t sep;
        // read element symbol (e.g. HE:NASA)
        if ((sep = line.find(':')) == string::npos)
            throw BasisSetFormatError(file, "':' not found", lineno);
        if (sep < 1 || sep > 2)
            throw BasisSetFormatError(file, "':' in wrong place", lineno);
        Element e = Element::getElement(line.substr(0, sep).c_str());

        // read comment line
        line = readLine(ifs, file, lineno);

        //read number of shells
        int nshell = readValue<int>(ifs, file, lineno);

        vector<ShellBasis> sb(nshell);

        // read L for each shell
        vector<int> L = readValues<int>(ifs, file, lineno, nshell);
        for (int i = 0;i < nshell;i++) sb[i].L = L[i];

        // read ncontr for each shell
        vector<int> ncontr = readValues<int>(ifs, file, lineno, nshell);
        for (int i = 0;i < nshell;i++) sb[i].ncontr = ncontr[i];

        // read nprim for each shell
        vector<int> nprim = readValues<int>(ifs, file, lineno, nshell);
        for (int i = 0;i < nshell;i++) sb[i].nprim = nprim[i];

        for (int i = 0;i < nshell;i++)
        {
            ShellBasis& b = sb[i];
            b.coefficients.resize(b.nprim*b.ncontr);

            b.exponents = readValues<double>(ifs, file, lineno, b.nprim);
            vector<double> coef = readValues<double>(ifs, file, lineno, b.nprim*b.ncontr);

            for (int j = 0;j < b.nprim;j++)
            {
                for (int k = 0;k < b.ncontr;k++)
                {
                    b.coefficients[j + k*b.nprim] = coef[k + j*b.ncontr];
                }
            }
        }

        atomBases[string(e.getName())] = sb;
    }
}

string BasisSet::readLine(istream& is, const string& file, int& lineno)
{
    string line;
    lineno++;
    if (!getline(is, line)) throw BasisSetFormatError(file, "premature EOF", lineno);
    return line;
}

void BasisSet::apply(Atom& atom, bool spherical, bool contaminants)
{
    string e(atom.getCenter().getElement().getName());
    vector<ShellBasis> v;
    vector<ShellBasis>::iterator it2;

		if (e == "Dummy" || e == "Ghost") return;

    map< string,vector<ShellBasis> >::iterator it = atomBases.find(e);
    if (it == atomBases.end())
    {
        throw BasisSetNotFoundError(e);
    }

    v = it->second;

    for (it2 = v.begin();it2 != v.end();++it2)
    {
        atom.addShell(Shell(atom.getCenter(), it2->L, it2->nprim, it2->ncontr,
                      spherical, contaminants, it2->exponents, it2->coefficients));
    }
}

void BasisSet::apply(Molecule& molecule, bool spherical, bool contaminants)
{
    vector<Atom>::iterator it;

    for (it = molecule.getAtomsBegin();it != molecule.getAtomsEnd();++it)
    {
        apply(*it, spherical, contaminants);
    }
}

}
}
