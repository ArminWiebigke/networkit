#include <deque>
#include <set>
#include <algorithm>
#include <cmath>

#include "DequeNumeric.h"

namespace deque_numeric {

bool compare(std::deque<double> &one, std::deque<double> &two) {
    if (one.size() != two.size())
        return false;

    for (int i = 0; i < int(one.size()); i++) {
        if (fabs(one[i] - two[i]) > 1e-7)
            return false;
    }
    return true;
}

double euclidean_norm(const std::deque<double> &a) {
    double norm = 0;
    for (double i : a)
        norm += i * i;
    return sqrt(norm);
}

int euclidean_normalize(std::deque<double> &a) {
    double norm = euclidean_norm(a);
    for (double &i : a)
        i /= norm;
    return 0;
}

double scalar_product(std::deque<double> &a, std::deque<double> &b) {
    double norm = 0;
    for (int i = 0; i < int(a.size()); i++)
        norm += a[i] * b[i];
    return norm;
}

int orthogonalize(std::deque<double> &a, std::deque<std::deque<double>> &M) {
    euclidean_normalize(a);
    for (auto &i : M) {
        double w = scalar_product(a, i);
        for (int j = 0; j < int(a.size()); j++)
            a[j] -= w * i[j];
    }
    euclidean_normalize(a);
    return 0;
}

int matrix_time_vector(std::deque<std::deque<double>> &Q, std::deque<double> &v,
                       std::deque<double> &new_s) {
    new_s.clear();
    for (auto &i : Q) {
        double n = 0;
        for (int j = 0; j < int(i.size()); j++)
            n += i[j] * v[j];

        new_s.push_back(n);
    }
    return 0;
}

void set_to_deque(const std::set<int> &s, std::deque<int> &a) {
    a.clear();
    for (int its : s)
        a.push_back(its);
}

void set_to_deque(const std::set<double> &s, std::deque<double> &a) {
    a.clear();
    for (double its : s)
        a.push_back(its);
}

void deque_to_set(const std::deque<double> &a, std::set<double> &s) {
    s.clear();
    for (double i : a)
        s.insert(i);
}

void deque_to_set(const std::deque<int> &a, std::set<int> &s) {
    s.clear();
    for (int i : a)
        s.insert(i);
}

void deque_to_set_app(const std::deque<int> &a, std::set<int> &s) {
    for (int i : a)
        s.insert(i);
}

double norm_one(const std::deque<double> &a) {
    double norm = 0;
    for (double i : a)
        norm += i;
    return norm;
}

int normalize_one(std::deque<double> &a) {
    double norm = norm_one(a);
    for (double &i : a)
        i /= norm;
    return 0;
}

double jaccard(std::set<int> &a1, std::set<int> &a2) {
    std::deque<int> group_intsec;
    std::deque<int> group_union;
    std::set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_intsec));
    std::set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_union));
    return double(group_intsec.size()) / group_union.size();
}

}
