#pragma once

#ifndef DEQUE_NUMERIC_HPP
#define DEQUE_NUMERIC_HPP

#include <deque>
#include <set>

namespace deque_numeric {

bool compare(std::deque<double> &one, std::deque<double> &two);

double euclidean_norm(const std::deque<double> &a);

int euclidean_normalize(std::deque<double> &a);

double scalar_product(std::deque<double> &a, std::deque<double> &b);

int orthogonalize(std::deque<double> &a, std::deque<std::deque<double>> &M);

int matrix_time_vector(std::deque<std::deque<double>> &Q, std::deque<double> &v,
                       std::deque<double> &new_s);

void set_to_deque(const std::set<int> &s, std::deque<int> &a);

void set_to_deque(const std::set<double> &s, std::deque<double> &a);

void deque_to_set(const std::deque<double> &a, std::set<double> &s);

void deque_to_set(const std::deque<int> &a, std::set<int> &s);

void deque_to_set_app(const std::deque<int> &a, std::set<int> &s);

double norm_one(const std::deque<double> &a);

int normalize_one(std::deque<double> &a);

double jaccard(std::set<int> &a1, std::set<int> &a2);

}

#endif //DEQUE_NUMERIC_HPP
