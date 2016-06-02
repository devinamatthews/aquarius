#ifndef _AQUARIUS_CC_CFOURGRAD_HPP_
#define _AQUARIUS_CC_CFOURGRAD_HPP_

#include "util/global.hpp"

#include "task/task.hpp"
#include "time/time.hpp"
#include "operator/2eoperator.hpp"
#include "input/molecule.hpp"

namespace aquarius
{
namespace cc
{

class CFOURGradient : public task::Task
{
    protected:
        input::Config config;

    public:
        CFOURGradient(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);

        void writeZMAT();

        void writeGENBAS();

        void execute(const Arena& arena, const string& module);

        void readIntegrals(const Arena& arena);

        void readIntegrals(ifstream& ifs, tensor::CTFTensor<double>& H, bool transpq, bool transrs);

        void readIntegrals(ifstream& ifs, tensor::CTFTensor<double>& H);

        void writeDensity();

        void writeDensity(ofstream& ofs, const tensor::CTFTensor<double>& D, bool transpq, bool transrs);

        void writeDensity(ofstream& ofs, const tensor::CTFTensor<double>& D);

        void clean();
};

}
}

#endif
