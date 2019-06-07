#include "Combinatorics.h"

double log_factorial(int num) {
    double log_result = 0;
    for (int i = 1; i <= num; i++)
        log_result += std::log(i);

    return (log_result);
}

double log_combination(int n, int k) {

    if (k == 0)
        return 0;

    if (n < k)
        return 0;

    if (n - k < k)
        k = n - k;

    double log_c = 0;
    for (int i = n - k + 1; i <= n; i++)
        log_c += log(i);

    for (int i = 1; i <= k; i++)
        log_c -= log(i);

    return log_c;
}

double binomial(int n, int x, double p) {
    //	returns the binomial distribution, n trials, x successes, p probability
    if (p == 0) {
        if (x == 0)
            return 1;
        else
            return 0;
    }

    if (p >= 1) {
        if (x == n)
            return 1;
        else
            return 0;
    }

    double log_b = 0;
    log_b += log_combination(n, x) + x * log(p) + (n - x) * log(1 - p);
    return (exp(log_b));
}

int binomial_cumulative(int n, double p, std::deque<double> &cum) {
    cum.clear();
    double c = 0;
    for (int i = 0; i <= n; i++) {
        c += binomial(n, i, p);
        cum.push_back(c);
    }
    return 0;
}

int powerlaw(int n, int min, double tau, std::deque<double> &cumulative) {
    cumulative.clear();
    double a = 0;
    for (double h = min; h < n + 1; h++)
        a += pow((1. / h), tau);
    double pf = 0;
    for (double i = min; i < n + 1; i++) {
        pf += 1 / a * pow((1. / (i)), tau);
        cumulative.push_back(pf);
    }
    return 0;
}

int distribution_from_cumulative(const std::deque<double> &cum, std::deque<double> &distr) {
    // cum is the cumulative, distr is set equal to the distribution
    distr.clear();
    double previous = 0;
    for (double i : cum) {
        distr.push_back(i - previous);
        previous = i;
    }
    return 0;
}

int cumulative_from_distribution(std::deque<double> &cum, const std::deque<double> &distr) {
    // cum is set equal to the cumulative, distr is the distribution
    cum.clear();
    double sum = 0;
    for (double i : distr) {
        sum += i;
        cum.push_back(sum);
    }
    return 0;
}

double poisson(int x, double mu) {
    return (exp(-mu + x * log(mu) - log_factorial(x)));
}

int shuffle_and_set(int *due, int dim) {
    // it sets due as a random sequence of integers from 0 to dim-1
    std::multimap<double, int> uno;
    for (int i = 0; i < dim; i++)
        uno.insert(std::make_pair(ran4(), i));

    std::multimap<double, int>::iterator it;

    int h = 0;
    for (it = uno.begin(); it != uno.end(); it++)
        due[h++] = it->second;

    return 0;
}

double compute_r(int x, int k, int kout, int m) {
    double r = 0;
    for (int i = x; i <= k; i++)
        r += binomial(k, i, double(kout) / double(m));
    return r;
}

double compute_hypergeometric(int i, int k, int kout, int m) {
    if (i > k || i > kout || k > m || kout > m)
        return 0;
    double prod = 1;
    std::deque<double> num;
    std::deque<double> den;

    if (add_factors(num, den, kout, i) == -1)
        return 0;

    if (add_factors(num, den, m - kout, k - i) == -1)
        return 0;

    if (add_factors(den, num, m, k) == -1)
        return 0;

    std::sort(num.begin(), num.end());
    std::sort(den.begin(), den.end());

    //prints(den);

    for (double h : den)
        if (h <= 0) {
            std::cerr << "denominator has zero or less (in the hypergeometric)" << std::endl;
            return 0;
        }

    for (double h : num)
        if (h <= 0) {
            std::cerr << "numerator has zero or less (in the hypergeometric)" << std::endl;
            return 0;
        }

    //cout<<"sizes: "<<num.size()<<" "<<den.size()<<endl;
    for (int idx = 0; idx < int(num.size()); idx++)
        prod = prod * num[idx] / den[idx];
    return prod;
}

int random_from_set(std::set<int> &s) {
    int pos1 = irand(s.size() - 1);
    auto it1 = s.begin();
    for (int i = 0; i < pos1; i++)
        it1++;
    return *it1;
}

int add_factors(std::deque<double> &num, std::deque<double> &den, int n, int k) {
    if (n < k)
        return -1;
    if (n - k < k)
        k = n - k;
    if (k == 0)
        return 0;
    for (int i = n - k + 1; i <= n; i++)
        num.push_back(double(i));
    for (int i = 1; i <= k; i++)
        den.push_back(double(i));
    return 0;
}
