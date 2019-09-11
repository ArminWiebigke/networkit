#include "Histograms.h"

int intlog_binning(std::deque<int> c, int number_of_bins, std::deque<double> &Xs,
                   std::deque<double> &Ys, std::deque<double> &var) {
    // this function is to make a log_histogram along the x axis and to compute y-errors
    Xs.clear();
    Ys.clear();
    var.clear();

    std::deque<int> d;
    for (int i : c)
        if (i > 0)
            d.push_back(i);

    c.clear();
    c = d;
    std::sort(c.begin(), c.end());
    int max = c[c.size() - 1];
    int min = c[0];
    std::deque<double> bins;
    double step = std::log(min);
    if (max == min)
        max++;

    double bin = (log(max) - log(min)) / number_of_bins;        // bin width
    while (step <= log(max) + 4 * bin) {
        bins.push_back(exp(step));
        step += bin;
    }

    std::deque<int> hist;
    std::deque<double> hist2;
    int index = 0;
    for (int i : c) {
        while (i - bins[index] > -1e-6) {
            index++;
            hist.push_back(0);
            hist2.push_back(0);
        }
        hist[index - 1]++;
        hist2[index - 1] += double(i);
    }

    std::deque<int> integers;
    index = 0;
    for (int i = min; i < bins[bins.size() - 1] - 1; i++) {
        while (i - bins[index] > -1e-6) {
            index++;
            integers.push_back(0);
        }
        integers[index - 1]++;
    }

    for (int i = 0; i < int(hist.size()); i++) {
        if (hist[i] > 0) {
            Xs.push_back(hist2[i] / hist[i]);
            double y = double(hist[i]) / (c.size() * integers[i]);
            Ys.push_back(y);
            var.push_back(
                    double(hist[i]) / (c.size() * c.size() * integers[i] * integers[i]));
        }
    }
    return 0;
}

void compute_quantiles(double q, std::deque<double> &y, std::deque<double> &qs) {
    int qv = cast_int((1 - q) / 2 * y.size());
    if (qv < 0)
        qv = 0;
    if (qv >= int(y.size()))
        qv = y.size() - 1;

    qs.push_back(y[qv]);

    qv = cast_int((1 + q) / 2 * y.size());
    if (qv < 0)
        qv = 0;
    if (qv >= int(y.size()))
        qv = y.size() - 1;

    qs.push_back(y[qv]);
}

int
log_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out, int number_of_bins) {        // c is the set od data, min is the lower bound, max is the upper one
    // c must be sorted and must be c[i]>0
    auto min = double(c[0]);
    auto max = double(c[0]);

    for (double i : c) {
        if (min > double(i))
            min = double(i);
        if (max < double(i))
            max = double(i);
    }

    std::deque<int> hist0;
    std::deque<double> hist;
    std::deque<double> hist2;
    std::deque<double> bins;
    double step = log(min);
    if (max == min)
        max++;

    double bin = (log(max) - log(min)) / number_of_bins;        // bin width

    while (step <= log(max) + 2 * bin) {
        bins.push_back(exp(step));
        hist0.push_back(0);
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
        hist0[index]++;
        hist[index] += c2[i];
        hist2[index] += double(c[i]) * c2[i];
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        double y = double(hist[i]) / (hist0[i]);
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

int histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out, int number_of_bins,
              double b1, double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    auto min = double(c[0]);
    auto max = double(c[0]);
    for (double i : c) {
        if (min > double(i))
            min = double(i);
        if (max < double(i))
            max = double(i);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;

    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist0;
    std::deque<double> hist;
    std::deque<double> hist2;

    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width
    while (step <= max + 2 * bin) {
        hist0.push_back(0);
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        auto data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist0[index]++;
            hist[index] += c2[i];
            hist2[index] += double(c[i] * c2[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        auto y = double(hist[i] / hist0[i] / bin);
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

int not_norm_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream &out,
                       int number_of_bins, double b1, double b2) {
    // this should be OK
    // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
    auto min = double(c[0]);
    auto max = double(c[0]);

    for (double i : c) {
        if (min > double(i))
            min = double(i);
        if (max < double(i))
            max = double(i);
    }

    min -= 1e-6;
    max += 1e-6;

    if (b1 != b2) {
        min = b1;
        max = b2;
    }
    if (max == min)
        max += 1e-3;

    std::deque<int> hist0;
    std::deque<double> hist;
    std::deque<double> hist2;
    double step = min;
    double bin = (max - min) / number_of_bins;        // bin width

    while (step <= max + 2 * bin) {
        hist0.push_back(0);
        hist.push_back(0);
        hist2.push_back(0);
        step += bin;
    }

    for (int i = 0; i < int(c.size()); i++) {
        auto data = double(c[i]);
        if (data > min && data <= max) {
            int index = int((data - min) / bin);
            hist0[index]++;
            hist[index] += c2[i];
            hist2[index] += double(c[i] * c2[i]);
        }
    }

    for (int i = 0; i < int(hist.size()) - 1; i++) {
        double x = hist2[i] / hist[i];
        auto y = double(hist[i] / hist0[i]);
        if (fabs(y) > 1e-10)
            out << x << "\t" << y << std::endl;
    }
    return 0;
}

void int_histogram(std::vector<int> &c, std::ostream &out) {
    std::map<int, double> hist;
    double freq = 1 / double(c.size());
    for (int & i : c) {
        auto itf = hist.find(i);
        if (itf == hist.end())
            hist.insert(std::make_pair(i, 1.));
        else
            itf->second++;
    }

    for (auto & it : hist)
        it.second = it.second * freq;
    prints(hist, out);
}

void int_histogram(std::deque<int> &c, std::ostream &out) {
    std::map<int, double> hist;
    double freq = 1 / double(c.size());
    for (int & i : c) {
        auto itf = hist.find(i);
        if (itf == hist.end())
            hist.insert(std::make_pair(i, 1.));
        else
            itf->second++;
    }

    for (auto & it : hist)
        it.second = it.second * freq;
    prints(hist, out);
}

void int_histogram(int c, std::map<int, int> &hist) {
    auto itf = hist.find(c);
    if (itf == hist.end())
        hist.insert(std::make_pair(c, 1));
    else
        itf->second++;
}

void int_histogram(int c, std::map<int, double> &hist, double w) {
    auto itf = hist.find(c);
    if (itf == hist.end())
        hist.insert(std::make_pair(c, w));
    else
        itf->second += w;
}

int print_cumulative(std::deque<double> &kws, const std::string &file, int number_of_step) {
    char buffer[100];
    cast_string_to_char(file, buffer);
    std::ofstream expout(buffer);
    sort(kws.begin(), kws.end());

    int step = (kws.size() - 1) / number_of_step;
    step = std::max(step, 1);

    for (int i = 0; i < int(kws.size()); i++)
        if (i % step == 0)
            expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
    return 0;
}

int print_cumulative(std::deque<int> &kws, const std::string &file, int number_of_step) {
    char buffer[100];
    cast_string_to_char(file, buffer);
    std::ofstream expout(buffer);
    sort(kws.begin(), kws.end());

    int step = (kws.size() - 1) / number_of_step;
    step = std::max(step, 1);

    for (int i = 0; i < int(kws.size()); i++)
        if (i % step == 0)
            expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
    return 0;
}

int print_cumulative(std::vector<double> &kws, const std::string &file, int number_of_step) {
    char buffer[100];
    cast_string_to_char(file, buffer);
    std::ofstream expout(buffer);
    sort(kws.begin(), kws.end());

    int step = (kws.size() - 1) / number_of_step;
    step = std::max(step, 1);

    for (int i = 0; i < int(kws.size()); i++)
        if (i % step == 0)
            expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
    return 0;
}

int print_cumulative(std::vector<int> &kws, const std::string &file, int number_of_step) {
    char buffer[100];
    cast_string_to_char(file, buffer);
    std::ofstream expout(buffer);
    sort(kws.begin(), kws.end());

    int step = (kws.size() - 1) / number_of_step;
    step = std::max(step, 1);
    for (int i = 0; i < int(kws.size()); i++)
        if (i % step == 0)
            expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
    return 0;
}

void int_histogram(const std::string &infile, const std::string &outfile) {
    // this makes a int_histogram of integers from a file
    char b[1000];
    cast_string_to_char(infile, b);

    std::ifstream ing(b);
    std::deque<int> H;
    int h;
    while (ing >> h)
        H.push_back(h);

    cast_string_to_char(outfile, b);
    std::ofstream outg(b);
    int_histogram(H, outg);
}

void int_histogram(const int &c, std::map<int, std::pair<int, double> > &hist, const int &w1,
                   const double &w2) {
    auto itf = hist.find(c);
    if (itf == hist.end())
        hist.insert(make_pair(c, std::make_pair(w1, w2)));
    else {
        itf->second.first += w1;
        itf->second.second += w2;
    }
}

void int_histogram(int c, std::map<int, std::pair<double, double> > &hist, double w1, double w2) {
    auto itf = hist.find(c);
    if (itf == hist.end())
        hist.insert(make_pair(c, std::make_pair(w1, w2)));
    else {
        itf->second.first += w1;
        itf->second.second += w2;
    }
}

void int_histogram(int c, std::map<int, std::pair<int, int> > &hist, int w1, int w2) {
    auto itf = hist.find(c);
    if (itf == hist.end())
        hist.insert(make_pair(c, std::make_pair(w1, w2)));
    else {

        itf->second.first += w1;
        itf->second.second += w2;
    }
}

/**
 *
 * @param hist Insert result here
 * @param c Neighbor node
 * @param w Weight of edge to neighbor
 */
void int_histogram(std::map<int, int> &hist, int c, int w) {
	if (hist.count(c))
		hist[c] += w;
	else
		hist[c] = w;
}
