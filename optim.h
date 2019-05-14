#ifndef _OPTIM_H_
#define _OPTIM_H_
#include <string>
#include <vector>
#include <armadillo>
#include <map>

#include "MCL.h"
#include "optim.h"

// needed by alglib
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

class lm_data {
    public:
        lm_data() { hatTheta_c.zeros(4,1);    }

    std::string         base;
    std::vector<double> errRateV;
    arma::mat           hatTheta_c;
    std::map<double, double> errRateFrac;
};

typedef std::pair<char, unsigned long> pCharUlong;

void function1_grad(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *opt_data);
void function_noGrad(const alglib::real_1d_array &x, double &func, void *opt_data);

#endif // _OPTIM_H_
