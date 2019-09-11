#pragma once

#ifndef HISTOGRAMS_HPP
#define HISTOGRAMS_HPP

#include <deque>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "Wsarray.h"
#include "Cast.h"
#include "Print.h"
#include "Combinatorics.h"

int intlog_binning(std::deque<int> c, int number_of_bins, std::deque<double> &Xs,
                   std::deque<double> &Ys,
                   std::deque<double> &var);

template<typename type>
int xybinning(std::deque<type> &c, std::deque<type> &d, int number_of_bins, std::deque<double> &xs,
              std::deque<double> &ys, std::deque<double> &var, std::deque<int> &nums) {
    // so, this function takes two datasets (c and d) and gathers the data in bin, takes xs and ys as the average in each bin, var is the variance of the y average
    // the difference with the same stuff called not_norm_histogram is that the other one averages x with y weights.
    xs.clear();
    ys.clear();
    var.clear();
    nums.clear();
    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }
    min -= 1e-6;
    max += 1e-6;

    if (max == min)
        max += 1e-3;

    std::deque<std::deque<double>> hist_x;        // x values in the bin
    std::deque<std::deque<double>> hist_y;        // y values in the bin

    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    std::deque<double> f;
    while (step <= max + 2 * bin) {
        hist_x.push_back(f);
        hist_y.push_back(f);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data >= min && data <= max) {
            int index = int((data - min) / bin);
            hist_x[index].push_back(double(c[i]));
            hist_y[index].push_back(double(d[i]));
        }
    }

    for (int i = 0; i < (int) hist_x.size() - 1; i++) {
        double x = average_func(hist_x[i]);
        double y = average_func(hist_y[i]);
        if (!hist_y[i].empty()) {
            xs.push_back(x);
            ys.push_back(y);
            var.push_back(variance_func(hist_y[i]) / double(hist_y[i].size()));
            nums.push_back(hist_y[i].size());
        }
    }

    for (double &i : var)
        if (i < 1e-8)
            i = 1e-8;

    return 0;
}

template<typename type>
int xybinning(std::deque<type> &c, std::deque<type> &d, int number_of_bins, std::deque<double> &xs,
              std::deque<double> &ys, std::deque<double> &var) {
    std::deque<int> nums;
    return xybinning(c, d, number_of_bins, xs, ys, var, nums);
}

void compute_quantiles(double q, std::deque<double> &y, std::deque<double> &qs);

template<typename type>
int xybinning_quantiles(std::deque<type> &c, std::deque<type> &d, int number_of_bins,
                        std::deque<double> &xs, std::deque<double> &ys, std::deque<double> &var,
                        std::deque<int> &nums, std::deque<std::deque<double>> &Mq, double qa,
                        double qb) {
    // so, this function takes two datasets (c and d) and gathers the data in bin, takes xs and ys as the average in each bin, var is the variance of the y average
    // the difference with the same stuff called not_norm_histogram is that the other one averages x with y weights.
    xs.clear();
    ys.clear();
    var.clear();
    nums.clear();
    Mq.clear();

    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;

    if (max == min)
        max += 1e-3;

    std::deque<std::deque<double>> hist_x;        // x values in the bin
    std::deque<std::deque<double>> hist_y;        // y values in the bin

    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    std::deque<double> f;
    while (step <= max + 2 * bin) {
        hist_x.push_back(f);
        hist_y.push_back(f);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data >= min && data <= max) {
            int index = int((data - min) / bin);
            hist_x[index].push_back(double(c[i]));
            hist_y[index].push_back(double(d[i]));
        }
    }

    for (int i = 0; i < (int) hist_x.size() - 1; i++) {
        double x = average_func(hist_x[i]);
        double y = average_func(hist_y[i]);
        if (!hist_y[i].empty()) {
            xs.push_back(x);
            ys.push_back(y);
            var.push_back(variance_func(hist_y[i]) / double(hist_y[i].size()));
            nums.push_back(hist_y[i].size());
            sort(hist_y[i].begin(), hist_y[i].end());
            std::deque<double> qs;
            compute_quantiles(qa, hist_y[i], qs);
            compute_quantiles(qb, hist_y[i], qs);
            Mq.push_back(qs);
        }
    }

    for (double &i : var)
        if (i < 1e-8)
            i = 1e-8;

    return 0;
}

template<typename type>
int log_histogram(std::deque<type> &c, std::ostream &out,
                  int number_of_bins) {        // c is the set od data, min is the lower bound, max is the upper one
    std::deque<type> d;
    for (int i = 0; i < int(c.size()); i++)
        if (c[i] > 0)
            d.push_back(c[i]);

    c.clear();
    c = d;

    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    std::deque<int> hist;
    std::deque<double> hist2;
    std::deque<double> bins;
    double step = log(min);
    if (max == min)
        max++;

    double bin = (log(max) - log(min)) / number_of_bins;        // bin width

    while (step <= log(max) + 2 * bin) {
        bins.push_back(exp(step));
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        int index = bins.size() - 1;
        for (int j = 0; j < int(bins.size()) - 1; j++)
            if ((fabs(double(c[i]) - bins[j]) < 1e-7) ||
                (double(c[i]) > bins[j] && double(c[i]) < bins[j + 1])) {
                // this could be done in a more efficient way
                index = j;
                break;
            }
        hist[index]++;
        hist2[index] += double(c[i]);
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double h1 = bins[i];
        double h2 = bins[i + 1];
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (c.size() * (h2 - h1));
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

int log_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out,
                  int number_of_bins);

template<typename type>
int histogram(std::vector<type> &c, std::ostream &out, int number_of_bins, double b1, double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist;
    std::deque<double> hist2;
    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width
    while (step <= max + 2 * bin) {
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist[index]++;
            hist2[index] += double(c[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (c.size() * bin);
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

template<typename type>
int not_norm_histogram_correlated(std::deque<type> &c, std::deque<type> &d, std::ostream &out,
                                  int number_of_bins, double b1, double b2) {
    // c is the x axis, d the y, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist;            // frequency in the bin
    std::deque<double> hist_x;        // x sum in the bin
    std::deque<double> hist_y;        // y sum in the bin
    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    while (step <= max + 2 * bin) {
        hist.push_back(0);
        hist_x.push_back(0);
        hist_y.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist[index]++;
            hist_x[index] += double(c[i]);
            hist_y[index] += double(d[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist_x[i] / hist[i];
        double y = hist_y[i] / hist[i];;
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

template<typename type>
int histogram(std::deque<type> &c, std::ostream &out, int number_of_bins, double b1, double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist;
    std::deque<double> hist2;

    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    while (step <= max + 2 * bin) {
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist[index]++;
            hist2[index] += double(c[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (c.size() * bin);
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

template<typename type>
int not_norm_histogram(std::vector<type> &c, std::ostream &out, int number_of_bins, double b1,
                       double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    double min = double(c[0]);
    double max = double(c[0]);

    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist;
    std::deque<double> hist2;
    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    while (step <= max + 2 * bin) {
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist[index]++;
            hist2[index] += double(c[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (c.size());
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

template<typename type>
int not_norm_histogram(std::deque<type> &c, std::ostream &out, int number_of_bins, double b1,
                       double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    double min = double(c[0]);
    double max = double(c[0]);
    for (int i = 0; i < int(c.size()); i++) {
        if (min > double(c[i]))
            min = double(c[i]);
        if (max < double(c[i]))
            max = double(c[i]);
    }

    min -= 1e-6;
    max += 1e-6;
    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist;
    std::deque<double> hist2;
    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    while (step <= max + 2 * bin) {
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        double data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist[index]++;
            hist2[index] += double(c[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (c.size());
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << "\t" << sqrt(hist[i]) / c.size() << std::endl;
    }
    return 0;
}

int histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out, int number_of_bins,
              double b1, double b2);

int
not_norm_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out,
                   int number_of_bins,
                   double b1, double b2);

void int_histogram(std::vector<int> &c, std::ostream &out);

void int_histogram(std::deque<int> &c, std::ostream &out);

void int_histogram(int c, std::map<int, int> &hist);

void int_histogram(int c, std::map<int, double> &hist, double w);

int print_cumulative(std::deque<double> &kws, const std::string &file, int number_of_step);

int print_cumulative(std::deque<int> &kws, const std::string &file, int number_of_step);

int print_cumulative(std::vector<double> &kws, const std::string &file, int number_of_step);

int print_cumulative(std::vector<int> &kws, const std::string &file, int number_of_step);

void int_histogram(const std::string &infile, const std::string &outfile);

void int_histogram(std::map<int, int> &hist, int c, int w);

void int_histogram(const int &c, std::map<int, std::pair<int, double> > &hist, const int &w1,
                   const double &w2);

void int_histogram(int c, std::map<int, std::pair<double, double> > &hist, double w1, double w2);

void int_histogram(int c, std::map<int, std::pair<int, int> > &hist, int w1, int w2);

#endif
