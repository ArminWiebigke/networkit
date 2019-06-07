#include "Wsarray.h"
#include <map>
#include <iostream>

Wsarray::Wsarray(int a) {
    position = 0;
    _size_ = a;
    l = new int[_size_];
    w = new std::pair<int, double>[_size_];
}

Wsarray::~Wsarray() {
    delete[] l;
    l = nullptr;
    delete[] w;
    w = nullptr;
}

std::pair<int, std::pair<int, double> > Wsarray::posweightof(int x) {
    int i = find(x);
    if (i == -1)
        return (std::make_pair(-1, std::make_pair(0, 0.)));
    return (std::make_pair(i, w[i]));
}

int Wsarray::find(int a) {
    int one = 0;
    int two = position - 1;

    if (position == 0)
        return -1;

    if (a < l[one] || a > l[two])
        return -1;

    if (a == l[one])
        return one;

    if (a == l[two])
        return two;

    while (two - one > 1) {
        int middle = (two - one) / 2 + one;
        if (a == l[middle])
            return middle;
        if (a > l[middle])
            one = middle;
        else
            two = middle;
    }
    return -1;
}

void Wsarray::push_back(int a, int bb, double b) {
    l[position] = a;
    w[position++] = std::make_pair(bb, b);
}

void Wsarray::freeze() {
    std::map<int, std::pair<int, double> > M;
    for (int i = 0; i < position; i++) {
        /*	//this is to sum up multiple entries
        map<int, double>::iterator itf=M.find(l[i]);
        if (itf==M.end())
            M.insert(std::make_pair(l[i], w[i]));
        else
            itf->second+=w[i];
        //*/
        M.insert(std::make_pair(l[i], std::make_pair(w[i].first, w[i].second)));
    }

    if (_size_ != int(M.size())) {
        delete[] l;
        l = nullptr;

        delete[] w;
        w = nullptr;

        _size_ = M.size();
        position = M.size();

        l = new int[_size_];
        w = new std::pair<int, double>[_size_];
    }

    int poi = 0;
    for (auto & itm : M) {

        l[poi] = itm.first;
        w[poi] = itm.second;
        poi++;
    }
}

void prints(Wsarray *a) {
    prints(a, std::cout);
}

void prints(Wsarray *a, std::ostream &pout) {
    for (int i = 0; i < a->size(); i++)
        std::cout << a->l[i] << "\t" << a->w[i].first << " " << a->w[i].second << std::endl;
    pout << std::endl;
}

void prints(Wsarray &a, std::ostream &pout) {
    for (int i = 0; i < a.size(); i++)
        std::cout << a.l[i] << "\t" << a.w[i].first << " " << a.w[i].second << std::endl;
    pout << std::endl;
}

void prints(Wsarray &a) {
    for (int i = 0; i < a.size(); i++)
        std::cout << a.l[i] << "\t" << a.w[i].first << " " << a.w[i].second << std::endl;
    std::cout << std::endl;
}
