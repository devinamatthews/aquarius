/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#include "basis.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;
using namespace aquarius::slide;

namespace aquarius
{
namespace input
{

BasisSet::ShellBasis::ShellBasis()
{
    L = 0;
    nprim = 0;
    ncontr = 0;
    exponents = NULL;
    coefficients = NULL;
}

BasisSet::ShellBasis::ShellBasis(const ShellBasis& other)
{
    L = other.L;
    nprim = other.nprim;
    ncontr = other.ncontr;
    exponents = new double[nprim];
    coefficients = new double[nprim*ncontr];
    copy(other.exponents, other.exponents+nprim, exponents);
    copy(other.coefficients, other.coefficients+nprim*ncontr, coefficients);
}

BasisSet::ShellBasis::~ShellBasis()
{
    delete[] exponents;
    delete[] coefficients;
}

BasisSet::ShellBasis& BasisSet::ShellBasis::operator=(ShellBasis other)
{
    swap(L, other.L);
    swap(nprim, other.nprim);
    swap(ncontr, other.ncontr);
    swap(exponents, other.exponents);
    swap(coefficients, other.coefficients);
    return *this;
}

BasisSet::BasisSet(const std::string& file)
{
    readBasisSet(file);
}

void BasisSet::readBasisSet(const string& file)
{
    ifstream ifs(file.c_str());
    string line;
    size_t sep;
    int lineno;
    int nshell;
    int i, j, k;
    istringstream iss;
    vector<ShellBasis> sb;
    ShellBasis* b;

    if (!ifs) throw BasisSetNotFoundError(file);

    for (lineno = 1;getline(ifs, line);lineno++)
    {
        // skip blank and comment lines
        if (line.find_first_not_of(" \t\n") == string::npos || line[0] == '!') continue;

        // read element symbol (e.g. HE:NASA)
        if ((sep = line.find(':')) == string::npos)
            throw BasisSetFormatError(file, "':' not found", lineno);
        if (sep < 1 || sep > 2)
            throw BasisSetFormatError(file, "':' in wrong place", lineno);
        Element e = Element::getElement(line.substr(0, sep).c_str());

        sb.clear();

        // read comment line
        line = readLine(ifs, file, lineno);

        //read number of shells
        nshell = readValue<int>(ifs, file, lineno);

        sb.resize(nshell);

        // read L for each shell
        int *L = new int[nshell];
        readValues<int>(ifs, file, lineno, nshell, L);
        for (i = 0;i < nshell;i++) sb[i].L = L[i];
        delete[] L;

        // read ncontr for each shell
        int *ncontr = new int[nshell];
        readValues<int>(ifs, file, lineno, nshell, ncontr);
        for (i = 0;i < nshell;i++) sb[i].ncontr = ncontr[i];
        delete[] ncontr;

        // read nprim for each shell
        int *nprim = new int[nshell];
        readValues<int>(ifs, file, lineno, nshell, nprim);
        for (i = 0;i < nshell;i++) sb[i].nprim = nprim[i];
        delete[] nprim;

        for (i = 0;i < nshell;i++)
        {
            b = &sb[i];
            b->exponents = new double[b->nprim];
            b->coefficients = new double[b->nprim*b->ncontr];
            double *coef = new double[b->nprim*b->ncontr];

            readValues<double>(ifs, file, lineno, b->nprim, b->exponents);
            readValues<double>(ifs, file, lineno, b->nprim*b->ncontr, coef);

            for (j = 0;j < b->nprim;j++)
            {
                for (k = 0;k < b->ncontr;k++)
                {
                    b->coefficients[j + k*b->nprim] = coef[k + j*b->ncontr];
                }
            }

            delete[] coef;
        }

        atomBases[string(e.getName())] = sb;
    }
}

template <typename T>
T BasisSet::readValue(istream& is, const string& file, int& lineno)
{
    T x;
    readValues<T>(is, file, lineno, 1, &x);
    return x;
}

template <typename T>
void BasisSet::readValues(istream& is, const string& file, int& lineno, int n, T* v)
{
    string line;
    istringstream iss;
    char c;

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
