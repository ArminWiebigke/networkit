#pragma once

#ifndef WSARRAY_HPP
#define WSARRAY_HPP

#include <iostream>

// this class is a container of int and double;
// it needs a preallocation when constructing the class;
// then you need to use the function push_back() to insert the numbers; this is the only way
// once you did, you can freeze the numbers; multiple entries are neglected; the preallocation can be reduced; you cannot use push_back anymore
// find returns the position of the integer you are looking for;
// posweightof returns a pair position-weight
// to access the integers use l[] and w[]
class Wsarray {

public:

    explicit Wsarray(int a);

    ~Wsarray();

    int find(int);

    std::pair<int, std::pair<int, double> > posweightof(int x);

    int size() { return position; };

    void freeze();

    void push_back(int, int, double);

    std::pair<int, double> *w;
    int *l;

private:

    int _size_;
    int position;
};

void prints(Wsarray &a);

void prints(Wsarray &a, std::ostream &pout);

void prints(Wsarray *a, std::ostream &pout);

void prints(Wsarray *a);

#endif



