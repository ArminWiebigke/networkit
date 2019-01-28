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

namespace NetworKit {

class UniformRandomSelector {
public:
    UniformRandomSelector();

    bool addElement();

    void reset();

private:
    count counter;
    std::mt19937 gen;
};

}

#endif //NETWORKIT_UNIFORMRANDOMSELECTOR_H
