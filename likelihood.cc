#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>

#include "main.h"
#include "likelihood.h"
#include <boost/math/distributions/chi_squared.hpp>

// needed by alglib
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "optimization.h"
#include "ap.h"

//#include "optim.h"
//using namespace arma;

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
    double l_c(0.0);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i];
        switch( base[i] ) {
            case 'A': l_c += log( (1-4*e/3) * theta(0) + e/3 ); break;
            case 'C': l_c += log( (1-4*e/3) * theta(1) + e/3 ); break;
            case 'G': l_c += log( (1-4*e/3) * theta(2) + e/3 ); break;
            case 'T': l_c += log( (1-4*e/3) * theta(3) + e/3 ); break;
            default: cerr << "unknown base in " << base << endl,exit(1);
        }
    }

    return l_c;
}

// composite score function: U_c(theta)
// return a column vector
//
alglib::real_1d_array composite_score (
        const string         &base,
        const vector<double> &errRateV,
        const alglib::real_1d_array  &theta )
{
   // mat U_c(4, 1, fill::zeros);
    alglib::real_1d_array U_c = "[0,0,0,0]";

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i];

        switch( base[i] ) {
            case 'A': U_c(0) -= (1-4*e/3) / ( (1-4*e/3)*theta(0) + e/3 );  break;
            case 'C': U_c(1) -= (1-4*e/3) / ( (1-4*e/3)*theta(1) + e/3 );  break;
            case 'G': U_c(2) -= (1-4*e/3) / ( (1-4*e/3)*theta(2) + e/3 );  break;
            case 'T': U_c(3) -= (1-4*e/3) / ( (1-4*e/3)*theta(3) + e/3 );  break;
            default: cerr << "unknown base in " << base << endl, exit(1);
        }
    }

    return U_c;
}

mCharDouble llh_genotype(const string &s, const string &q, const Option &opt)
{
    mCharDouble ntP; // nt => pvalue
    boost::math::chi_squared X2_dist(1);

    mCharUlong fr;
    double depth( s.size() );

    string new_s(""), new_q("");
    for ( size_t i(0); i != depth; i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' )  continue;
        fr[ s[i] ]++;
        new_s += s[i];
        new_q += q[i];
    }

    if ( fr.empty() )  return ntP;

    vector<double> errV = quaToErrorRate(new_q, opt);

    fn_data data;

    data.base = new_s;
    data.errRateV = errV;

    // four allele maximize
    alglib::real_1d_array alg_x = "[0.25,0.25,0.25,0.25]";
    alglib::real_2d_array c = "[[1,1,1,1,1]]";
    alglib::integer_1d_array ct = "[0]";
    alglib::minbleicstate state;
    alglib::minbleicreport rep;
    alglib::real_1d_array bndl = "[0,0,0,0]";
    alglib::real_1d_array bndu = "[1,1,1,1]";

    double epsg = 0.000001;
    double epsf = 0;
    double epsx = 0;
    alglib::ae_int_t maxits = 0;

    alglib::minbleiccreate(alg_x, state);
    alglib::minbleicsetlc(state, c, ct);
    alglib::minbleicsetbc(state, bndl, bndu);
    alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
    alglib::minbleicoptimize(state, function1_grad, NULL, &data );
    alglib::minbleicresults(state, alg_x, rep);

    //printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    if ( opt.debug ) printf("%s\n", alg_x.tostring(20).c_str());

    double cl_4 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
    if ( opt.debug ) cout << "l_c_hatTheta_c: " << setprecision(20) << cl_4 << endl;
/*
    alglib::real_1d_array fn_x = "[0.25,0.25,0.25,0.25]";

    double lm_4(0.0);
    try {
        double diffstep = 1.0e-6;
        alglib::minbleiccreatef(fn_x, diffstep, state);
        alglib::minbleicsetlc(state, c, ct);
        alglib::minbleicsetbc(state, bndl, bndu);
        alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
        alglib::minbleicoptimize(state, function_noGrad, NULL, &d );
        alglib::minbleicresults(state, fn_x, rep);
    }
    catch (alglib::ap_error & e) {
        lm_4 = cl_4;
 //       cerr << "catch exception & throw " << e.msg << " at: " << data.base << endl;
    }

    //printf("%d\n", int(rep.terminationtype));
    //printf("%s\n", fn_x.tostring(2).c_str());

    arma::mat fnt(4,1,arma::fill::zeros);
    for ( size_t i(0); i<4; i++ ) fnt(i) =  fn_x[i];

    if ( lm_4 == 0.0)
        lm_4 = get_l_M_theta_origin( d.base, d.errRateV, fnt, d.hatTheta_c, d.errRateFrac );

    if ( int(rep.terminationtype) != 4 ) lm_4 = cl_4;

    //cout << "l_M_theta: " << lm_4 << endl;
*/
//    cout << "~~~ test alglib at frac(T) = 0 ~~~" << endl;

    map<char, string> init_V;
    map<char, string> bndu_V;


    init_V['A'] = "[0.0,0.25,0.25,0.5]";
    init_V['C'] = "[0.25,0.0,0.25,0.5]";
    init_V['G'] = "[0.25,0.25,0.0,0.5]";
    init_V['T'] = "[0.25,0.25,0.5,0.0]";

    bndu_V['A'] = "[0,1,1,1]"; //"[1e-9,1,1,1]";
    bndu_V['C'] = "[1,0,1,1]"; //"[1,1e-9,1,1]";
    bndu_V['G'] = "[1,1,0,1]"; //"[1,1,1e-9,1]";
    bndu_V['T'] = "[1,1,1,0]"; //"[1,1,1,1e-9]";

    for ( mCharUlong::const_iterator it = fr.begin(); it != fr.end(); it++ )
    {
        alglib::real_1d_array alg_x = init_V[ it->first ].c_str();
        alglib::real_2d_array c = "[[1,1,1,1,1]]";
        alglib::integer_1d_array ct = "[0]";
        alglib::minbleicstate state;
        alglib::minbleicreport rep;
        alglib::real_1d_array bndl = "[0.0, 0.0, 0.0, 0.0]";
        alglib::real_1d_array bndu = bndu_V[ it->first ].c_str();


        double epsg = 0.000001;
        double epsf = 0;
        double epsx = 0;
        alglib::ae_int_t maxits = 0;

        alglib::minbleiccreate(alg_x, state);

        alglib::minbleicsetlc(state, c, ct);
        alglib::minbleicsetbc(state, bndl, bndu);
        alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);
        alglib::minbleicoptimize(state, function1_grad, NULL, &data );
        alglib::minbleicresults(state, alg_x, rep);

        //printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
        if ( opt.debug ) printf("%s\n", alg_x.tostring(20).c_str());

        double cl_3 = composite_LogLikelihood( data.base, data.errRateV, alg_x );
        if ( opt.debug ) cout << "3 l_c_hatTheta_c: " << cl_3 << endl;
/*
        alglib::real_1d_array fn_x = init_V[ it->first ].c_str();

        double lm_3(0.0);
        try {
            double diffstep = 1.0e-6;
            alglib::minbleiccreatef(fn_x, diffstep, state);
            alglib::minbleicsetlc(state, c, ct);
            alglib::minbleicsetbc(state, bndl, bndu);
            alglib::minbleicsetcond(state, epsg, epsf, epsx, maxits);

            alglib::minbleicoptimize(state, function_noGrad, NULL, &d );

            alglib::minbleicresults(state, fn_x, rep);
        }
        catch (alglib::ap_error & e) {
            lm_3 = cl_3;
//            cerr << "catch exception & throw " << e.msg << endl;
        }

        //printf("%d\n", int(rep.terminationtype));
        //printf("%s\n", fn_x.tostring(2).c_str());
        //cout << fn_x[3] << endl;

        arma::mat fnt(4,1,arma::fill::zeros);
        for ( size_t i(0); i<4; i++ ) fnt(i) =  fn_x[i];

        if (lm_3 == 0.0)
            lm_3 = get_l_M_theta_origin( d.base, d.errRateV, fnt, d.hatTheta_c, d.errRateFrac );

        if ( int(rep.terminationtype) != 4 ) lm_3 = cl_3;

        //cout << "3 l_M_theta: " << lm_3 << endl;
*/
//        boost::math::chi_squared X2_dist(1);

//        cout << "cl: " << 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3) ) << endl;
  //      cout << "lm: " << 1 - boost::math::cdf(X2_dist, 2*(lm_4 - lm_3) ) << endl;
        if ( cl_4 - cl_3 > opt.lhrGapCutoff )
            ntP[ it->first ] = 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3) );
    //
  //      if ( lm_4 - lm_3 > opt.lhrGapCutoff )
    //        ntP[ it->first ] = 1 - boost::math::cdf(X2_dist, 2*(lm_4 - lm_3) );
    }

    return ntP;

}

void function1_grad(const alglib::real_1d_array &x, double &func, alglib::real_1d_array &grad, void *opt_data)
{
    fn_data* objfn_data = reinterpret_cast<fn_data*>(opt_data);
    const std::string &base = objfn_data->base;
    const std::vector<double> &errRateV = objfn_data->errRateV;

    func = -composite_LogLikelihood( base, errRateV, x);
    grad = composite_score( base, errRateV, x );
}

