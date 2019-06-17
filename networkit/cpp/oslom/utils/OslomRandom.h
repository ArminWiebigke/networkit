#pragma once

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <deque>
#include <set>

double ran2(long *idum);

double ran4(bool t, long s);

/**
 * @return a random value in the interval [0, 1]
 */
double ran4();

void srand4();

void srand5(int rank);

int irand(int n);

void srand_file();

int configuration_model(std::deque <std::set<int>> &en, std::deque<int> &degrees);

#endif //RANDOM_HPP
