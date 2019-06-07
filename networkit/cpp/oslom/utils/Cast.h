#pragma once

#ifndef CAST_HPP
#define CAST_HPP

#include <string>
#include <deque>

bool cast_string_to_double(std::string &b, double &h);

double cast_string_to_double(std::string &b);

inline int cast_int(double u) {
    return std::lround(u);
}

int cast_string_to_char(const std::string &file_name, char *b);

bool cast_string_to_doubles(std::string &b, std::deque<double> &v);

bool cast_string_to_doubles(std::string &b, std::deque<int> &v);

bool separate_strings(std::string &b, std::deque <std::string> &v);

double approx(double a, int digits);

#endif //CAST_HPP
