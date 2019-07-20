//
// Created by armin on 7/20/19.
//

#ifndef NETWORKIT_PARSESTRING_HPP
#define NETWORKIT_PARSESTRING_HPP

#include <string>

namespace Aux {

/**
 * Convert a string to a floating point number.
 * @param str
 * @param decimalChar The character representing the decimal point, either ',' or '.'
 * @return
 */
double stringToDouble(std::string str, char decimalChar = '.');

}

#endif //NETWORKIT_PARSESTRING_HPP
