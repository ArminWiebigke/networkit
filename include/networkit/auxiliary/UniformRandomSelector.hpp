/*
 * UniformRandomSelector.hpp
 *
 * Created: 2019-01-21
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_UNIFORMRANDOMSELECTOR_HPP
#define NETWORKIT_UNIFORMRANDOMSELECTOR_HPP

#include <random>

#include <networkit/Globals.hpp>

namespace Aux {

class UniformRandomSelector {
public:
    UniformRandomSelector();

    bool addElement();

    void reset();

private:
    NetworKit::count counter;
};

}

#endif //NETWORKIT_UNIFORMRANDOMSELECTOR_HPP
