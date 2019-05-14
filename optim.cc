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

void function1_grad(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *opt_data)
{
    fn_data* objfn_data = reinterpret_cast<fn_data*>(opt_data);
    const std::string &base = objfn_data->base;
    const std::vector<double> &errRateV = objfn_data->errRateV;

    arma::mat theta(4,1,arma::fill::zeros);
    for ( size_t i(0); i < 4; i++  )  theta(i) = x[i];

    func = -composite_LogLikelihood( base, errRateV, theta);
    arma::mat score = -composite_score( base, errRateV, theta );

    for ( size_t i(0); i < 4; i++ )   grad[i] = score(i);
}

void function_noGrad(const alglib::real_1d_array &x, double &func, void *opt_data)
{
    lm_data* dat= reinterpret_cast<lm_data*>(opt_data);

    arma::mat theta(4,1,arma::fill::zeros);
    for ( size_t i(0); i < 4; i++  )  theta(i) = x[i];

    func = -get_l_M_theta_origin( dat->base, dat->errRateV, theta, dat->hatTheta_c, dat->errRateFrac );
}
