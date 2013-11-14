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

#ifndef _AQUARIUS_INPUT_BASIS_HPP_
#define _AQUARIUS_INPUT_BASIS_HPP_

#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "integrals/shell.hpp"
#include "molecule.hpp"

namespace aquarius
{
namespace input
{

class BasisSetNotFoundError : public std::runtime_error
{
    public:
        BasisSetNotFoundError(const std::string& element) : runtime_error(element) {}
};

class BasisSetFormatError : public std::runtime_error
{
    private:
        static std::string buildString(const std::string& file, const std::string& what_arg, const int lineno)
        {
            std::ostringstream os;
            os << file << ": " << what_arg << ": line " << lineno;
            return os.str();
        }

    public:
        BasisSetFormatError(const std::string& file, const std::string& what_arg, const int lineno)
        : runtime_error(buildString(file, what_arg, lineno)) {}
};

class BasisSet
{
    private:
        struct ShellBasis
        {
            std::vector<double> exponents;
            std::vector<double> coefficients;
            int nprim;
            int ncontr;
            int L;
        };

        std::map< std::string,std::vector<ShellBasis> > atomBases;

        void readBasisSet(const std::string& file);

        template <typename T>
        T readValue(std::istream& is, const std::string& file, int& lineno)
        {
            std::string line;
            std::istringstream iss;
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
        std::vector<T> readValues(std::istream& is, const std::string& file, int& lineno, int n)
        {
            std::string line;
            std::istringstream iss;
            char c;
            std::vector<T> v(n);

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

        std::string readLine(std::istream& is, const std::string& file, int& lineno);

    public:
        BasisSet() {}

        BasisSet(const std::string& file);

        void apply(Atom& atom, bool spherical = true, bool contaminants = false);

        void apply(Molecule& molecule, bool spherical = true, bool contaminants = false);
};

}
}

#endif
