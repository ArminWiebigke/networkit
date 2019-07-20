//
// Created by armin on 7/20/19.
//
#include <iostream>
#include <algorithm>

#include "ParseString.h"

namespace Aux {

double stringToDouble(std::string str, char decimalChar) {
    std::string old = str;
    if (decimalChar == '.') {
        if (std::stod("0.3") == 0) {
            std::replace(str.begin(), str.end(), ',', '#');
            std::replace(str.begin(), str.end(), '.', ',');
            std::replace(str.begin(), str.end(), '#', '.');
        }
    } else if (decimalChar == ',') {
        if (std::stod("0,3") == 0) {
            std::replace(str.begin(), str.end(), '.', '#');
            std::replace(str.begin(), str.end(), ',', '.');
            std::replace(str.begin(), str.end(), '#', ',');
        }
    }
    return std::stod(str);
}

}