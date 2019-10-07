/*
 * ParseString.cpp
 *
 * Created: 2019-07-20
 * Author: Armin Wiebigke
 */
#include <algorithm>

#include "ParseString.h"

namespace Aux {

void swapChars(std::string &str, char char1, char char2, char tmpChar);

double stringToDouble(std::string str, char decimalChar) {
	std::string old = str;
	// Use decimal separator of current locale
	if ((decimalChar == '.' && std::stod("0.3") == 0) ||
	    (decimalChar == ',' && std::stod("0,3") == 0)) {
		swapChars(str, '.', ',', '#');
	}
	return std::stod(str);
}

void swapChars(std::string &str, char char1, char char2, char tmpChar) {
	std::replace(str.begin(), str.end(), char1, tmpChar);
	std::replace(str.begin(), str.end(), char2, char1);
	std::replace(str.begin(), str.end(), tmpChar, char2);
}

}