#pragma once

#ifndef PRINT_HPP
#define PRINT_HPP

#include <iostream>
#include <cmath>
#include <map>
#include <deque>
#include <vector>
#include <fstream>

#include "Cast.h"

void cherr();

void cherr(double a);

void cherr(double a, double ee);

template<typename uno, typename due>
void prints(std::pair<uno, due> &sq, std::ostream &out) {
    out << sq.first << "\t" << sq.second << std::endl;
}

template<typename uno, typename due>
void prints(std::pair<uno, due> &sq) {
    std::cout << sq.first << "\t" << sq.second << std::endl;
}

template<typename uno, typename due>
void prints(std::map<uno, due> &sq, std::ostream &out) {
    typename std::map<uno, due>::iterator it = sq.begin();
    while (it != sq.end()) {
        out << it->first << "\t" << it->second << std::endl;
        it++;
    }
    out << std::endl;
}

template<typename uno, typename due>
void prints(std::multimap<uno, due> &sq, std::ostream &out) {
    typename std::map<uno, due>::iterator it = sq.begin();
    while (it != sq.end()) {
        out << it->first << "\t" << it->second << std::endl;
        it++;
    }
    out << std::endl;
}

template<typename Seq>
void prints(Seq &sq, std::ostream &out) {
    typename Seq::iterator it = sq.begin();
    while (it != sq.end())
        out << *(it++) << "\t";
    out << std::endl;
}

template<typename type_>
void prints(type_ *a, int b) {
    for (int i = 0; i < b; i++)
        std::cout << a[i] << " ";
    std::cout << std::endl;
}

template<typename T, template<typename> class C>
void printm(C<T> &c, std::ostream &out) {
    typename C<T>::iterator it = c.begin();
    while (it != c.end()) {
        prints(*it, out);
        it++;
    }
    out << std::endl;
}

template<typename uno, typename due>
void prints(std::map<uno, due> &sq) {
    typename std::map<uno, due>::iterator it = sq.begin();
    while (it != sq.end()) {
        std::cout << it->first << "\t" << it->second << std::endl;
        it++;
    }
    std::cout << std::endl;
}

template<typename uno, typename due>
void prints(std::multimap<uno, due> &sq) {
    typename std::map<uno, due>::iterator it = sq.begin();
    while (it != sq.end()) {
        std::cout << it->first << "\t" << it->second << std::endl;
        it++;
    }
    std::cout << std::endl;
}

template<typename uno, typename due>
void prints(std::deque <uno> &a, std::deque <due> &b) {
    for (int i = 0; i < a.size(); i++)
        std::cout << a[i] << " " << b[i] << std::endl;
}

template<typename uno, typename due>
void prints(std::deque <uno> &a, std::deque <due> &b, std::ostream &out) {
    for (int i = 0; i < a.size(); i++)
        out << a[i] << " " << b[i] << std::endl;
}

template<typename Seq>
void prints(Seq &sq) {
    typename Seq::iterator it = sq.begin();
    while (it != sq.end())
        std::cout << *(it++) << "\t";
    std::cout << std::endl;
}

template<typename type>
void prints(const std::deque <type> &sq) {
    for (int i = 0; i < sq.size(); i++)
        std::cout << sq[i] << "\t";
    std::cout << std::endl;
}

template<typename type>
void prints(const std::vector <type> &sq) {
    for (int i = 0; i < sq.size(); i++)
        std::cout << sq[i] << "\t";
    std::cout << std::endl;
}

template<typename type>
void printm(std::deque <type> &M) {
    for (int i = 0; i < M.size(); i++)
        prints(M[i]);
    std::cout << std::endl;
}

template<typename type>
void printm(std::vector <type> &M) {
    for (int i = 0; i < M.size(); i++)
        prints(M[i]);
    std::cout << std::endl;
}

template<typename type>
void get_data_from_file(const std::string& s, std::deque <type> &a1, int col) {
    // default will be col=1
    char b[1000];
    cast_string_to_char(s, b);
    std::ifstream lin(b);

    a1.clear();
    col--;

    std::string sas;
    while (getline(lin, sas)) {
        std::deque<double> doubles;
        cast_string_to_doubles(sas, doubles);
        if ((int) doubles.size() > col) {
            a1.push_back(doubles[col]);
        }
    }
}

template<typename type>
void get_data_from_file(std::string s, std::deque <type> &a1) {
    get_data_from_file(s, a1, 1);
}

template<typename type>
void
get_data_from_file(const std::string& s, std::deque <type> &a1, std::deque <type> &a2, int col1,
        int col2) {
    // default will be col1=1, col2=2
    char b[1000];
    cast_string_to_char(s, b);
    std::ifstream lin(b);

    a1.clear();
    a2.clear();
    col1--;
    col2--;

    std::string sas;
    while (getline(lin, sas)) {
        std::deque<double> doubles;
        cast_string_to_doubles(sas, doubles);
        if ((int) doubles.size() > col2) {
            a1.push_back(doubles[col1]);
            a2.push_back(doubles[col2]);
        }
    }
}

template<typename type>
void get_data_from_file(std::string s, std::deque <type> &a1, std::deque <type> &a2) {
    get_data_from_file(s, a1, a2, 1, 2);
}

void get_data_from_file_string(const std::string& s, std::deque <std::string> &a1, int col);

#endif
