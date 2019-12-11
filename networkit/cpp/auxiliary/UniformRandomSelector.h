/*
 * UniformRandomSelector.h
 *
 * Created: 2019-01-21
 * Author: Armin Wiebigke
 */

#ifndef NETWORKIT_UNIFORMRANDOMSELECTOR_H
#define NETWORKIT_UNIFORMRANDOMSELECTOR_H

#include <random>

#include "../Globals.h"

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

#endif //NETWORKIT_UNIFORMRANDOMSELECTOR_H
