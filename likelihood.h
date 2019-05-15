#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "ap.h"

struct fn_data {
    std::string         base;
    std::vector<double> errRateV;
};

//typedef std::pair<char, unsigned long> pCharUlong;

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta );

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
alglib::real_1d_array composite_score (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta );

mCharDouble llh_genotype(const string &s, const string &q, const Option &opt);

void function1_grad(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *opt_data);

#endif // LIKELIHOOD_H_
