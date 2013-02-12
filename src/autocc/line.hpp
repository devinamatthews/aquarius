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
 * ARE DISCLAIMED. IN NO EVENT SHALL EDGAR SOLOMONIK BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_AUTOCC_LINE_HPP_
#define _AQUARIUS_AUTOCC_LINE_HPP_

#include <ostream>
#include <string>
#include <vector>
#include <stdint.h>

namespace aquarius
{
namespace autocc
{

class Manifold;
class Line;

std::ostream& operator<<(std::ostream& out, const std::vector<Line>& v);
std::ostream& operator<<(std::ostream& out, const Line& l);
std::ostream& operator<<(std::ostream& out, const Manifold& m);

class Manifold
{
    private:
        void parse(const std::string& s);

    public:
        int np, nh;
        static const Manifold MIN_VALUE;
        static const Manifold MAX_VALUE;

        Manifold() : np(0), nh(0) {}

        Manifold(const Manifold& other) : np(other.np), nh(other.nh) {}

        Manifold(const char *s);

        Manifold(const std::string& s);

        Manifold(int ex) : np(ex), nh(ex) {}

        Manifold(int np, int nh) : np(np), nh(nh) {}

        Manifold& operator+=(const Manifold& other);

        Manifold& operator-=(const Manifold& other);

        Manifold operator+(const Manifold& other) const;

        Manifold operator-(const Manifold& other) const;

        bool operator<(const Manifold& other) const;

        bool operator==(const Manifold& other) const;
};

enum
{
    HOLE     = 0x04,
    PARTICLE = 0x08,
    INTERNAL = 0x01,
    EXTERNAL = 0x02,
    ALPHA    = 0x10,
    BETA     = 0x20
};

class Line
{
    friend std::ostream& operator<<(std::ostream& out, const Line& l);
    friend std::ostream& operator<<(std::ostream& out, const std::vector<Line>& v);

    protected:
        uint16_t index;
        static const unsigned int INDEX_MASK = 0x1FFF;
        static const unsigned int TYPE_MASK = 0xE000;
        static const unsigned int INDEX_MAX = 0x1FFF;
        enum {_HOLE = 0x4000, _INTERNAL = 0x2000, _ALPHA = 0x8000};

        Line(uint16_t index) : index(index) {}

    public:
        Line() : index(0) {}

        Line(const std::string& label);

        Line(char label);

        Line(unsigned int index, unsigned int type);

        static std::vector<Line> parse(std::string s);

        //operator unsigned int() const;

        Line& operator=(const Line& line);

        //Line& operator=(unsigned int index);

        bool operator<(const Line& other) const;

        bool operator==(const Line& other) const;

        int getIndex() const;

        int toInt() const;

        bool isExternal() const;

        bool isInternal() const;

        bool isHole() const;

        bool isParticle() const;

        bool isAlpha() const;

        bool isBeta() const;

        Line toAlpha() const;

        Line toBeta() const;
};

struct isExternal : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isExternal(); }
};

struct isInternal : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isInternal(); }
};

struct isHole : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isHole(); }
};

struct isParticle : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isParticle(); }
};

struct isAlpha : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isAlpha(); }
};

struct isBeta : std::unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isBeta(); }
};

template<int type> struct isType : std::unary_function<Line,bool> {};

template<> struct isType<HOLE> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole(); }
};

template<> struct isType<PARTICLE> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle(); }
};

template<> struct isType<INTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isInternal(); }
};

template<> struct isType<EXTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isExternal(); }
};

template<> struct isType<ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isAlpha(); }
};

template<> struct isType<BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isBeta(); }
};

template<> struct isType<HOLE+INTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isInternal(); }
};

template<> struct isType<HOLE+EXTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isExternal(); }
};

template<> struct isType<PARTICLE+INTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isInternal(); }
};

template<> struct isType<PARTICLE+EXTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isExternal(); }
};

template<> struct isType<HOLE+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isAlpha(); }
};

template<> struct isType<HOLE+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isBeta(); }
};

template<> struct isType<PARTICLE+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isAlpha(); }
};

template<> struct isType<PARTICLE+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isBeta(); }
};

template<> struct isType<ALPHA+INTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isAlpha() && line.isInternal(); }
};

template<> struct isType<ALPHA+EXTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isAlpha() && line.isExternal(); }
};

template<> struct isType<BETA+INTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isBeta() && line.isInternal(); }
};

template<> struct isType<BETA+EXTERNAL> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isBeta() && line.isExternal(); }
};

template<> struct isType<HOLE+INTERNAL+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isInternal() && line.isAlpha(); }
};

template<> struct isType<HOLE+INTERNAL+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isInternal() && line.isBeta(); }
};

template<> struct isType<HOLE+EXTERNAL+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isExternal() && line.isAlpha(); }
};

template<> struct isType<HOLE+EXTERNAL+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isHole() && line.isExternal() && line.isBeta(); }
};

template<> struct isType<PARTICLE+INTERNAL+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isInternal() && line.isAlpha(); }
};

template<> struct isType<PARTICLE+INTERNAL+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isInternal() && line.isBeta(); }
};

template<> struct isType<PARTICLE+EXTERNAL+ALPHA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isExternal() && line.isAlpha(); }
};

template<> struct isType<PARTICLE+EXTERNAL+BETA> : std::unary_function<Line,bool>
{
        bool operator()(const Line& line) const { return line.isParticle() && line.isExternal() && line.isBeta(); }
};

}
}

#endif
