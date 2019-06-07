#pragma once

#ifndef DATA_TYPES_HPP
#define DATA_TYPES_HPP

#include <map>

typedef std::multimap<double, std::pair<int, double> > CupDataStruct; // Fitness/Score, (node, (score)interval)

#endif //DATA_TYPES_HPP
