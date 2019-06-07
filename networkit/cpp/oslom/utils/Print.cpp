#include "Print.h"

void cherr() {
    std::cerr << "the check failed" << std::endl;
    int e;
    std::cin >> e;
}

void cherr(double a) {
    std::cerr << "the check failed because of " << a << std::endl;
    int e;
    std::cin >> e;
}

void cherr(double a, double ee) {
    if (fabs(a) > ee) {
        std::cerr << "the check failed because of " << a << std::endl;
        int e;
        std::cin >> e;
    }
}

void get_data_from_file_string(const std::string &s, std::deque<std::string> &a1, int col) {
    // default will be col=1
    char b[1000];
    cast_string_to_char(s, b);
    std::ifstream lin(b);

    a1.clear();
    col--;

    std::string sas;
    while (getline(lin, sas)) {
        std::deque <std::string> v;
        separate_strings(sas, v);
        //prints(v);
        if (int(s.size()) > col) {
            a1.push_back(v[col]);
        }
    }
}
