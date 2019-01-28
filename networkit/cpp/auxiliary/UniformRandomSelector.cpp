/*
 * UniformRandomSelector.cpp
 *
 * Created: 2019-01-21
 * Author: Armin Wiebigke
 */

#include <random>
#include <ctime>

#include "UniformRandomSelector.h"

namespace NetworKit {

UniformRandomSelector::UniformRandomSelector() : counter(1) {
    gen.seed(3249029);
}

bool UniformRandomSelector::addElement() {
    ++counter;
    std::uniform_int_distribution<count> dist(0, counter - 1);
    return (dist(gen) == 0);
}

void UniformRandomSelector::reset() {
    counter = 1;
}

}