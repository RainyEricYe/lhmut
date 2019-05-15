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

string initAlleleFreq (
        mCharUlong       &fr,
        const double     &depth,
        const char       &except_b )
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        b == except_b ? s << 0.0 : s << fr[b] / (depth - fr[except_b]);
        b == 'T' ? s << ']' : s << ',';
        if ( b== 'T' ) break;
    }
    return s.str();
}

string _upBoundary(const char except_b)
{
    ostringstream s;
    s << '[';

    for ( auto b : "ACGT" ) {
        b == except_b ? s << 0 : s << 1;
        b == 'T' ? s << ']' : s << ',';
        if ( b == 'T' ) break;
    }
    return s.str();
}


mCharDouble llh_genotype(const string &s, const string &q, const Option &opt)
{
    mCharDouble ntP; // nt => pvalue
    boost::math::chi_squared X2_dist(1);

    mCharUlong fr;
    double depth(0.0);
    for ( auto b : "ACGTN" ) {
        fr[b] = 0;
        if ( b == 'N' ) break;
    }
    // fr['N'] should always be 0

    string new_s(""), new_q("");
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' )  continue;
        fr[ s[i] ]++;
        new_s += s[i];
        new_q += q[i];
        depth++;
    }

    vector<double> errV = quaToErrorRate(new_q, opt);

    fn_data data;
    data.base = new_s;
    data.errRateV = errV;

    // four allele maximize
    //alglib::real_1d_array alg_x = "[0.25,0.25,0.25,0.25]";
    string AFstr = initAlleleFreq(fr, depth, 'N');
    alglib::real_1d_array alg_x = AFstr.c_str();
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

    map<char, string> init_V;
    map<char, string> bndu_V;

    for ( auto b : "ACGT" ) {
        init_V[b] = initAlleleFreq(fr, depth, b);
        bndu_V[b] = _upBoundary(b);
        if ( b == 'T' ) break;
    }

    for ( mCharUlong::const_iterator it = fr.begin(); it != fr.end(); it++ )
    {
        if ( it->second < opt.minSupOnEachStrand || it->second/depth < opt.minFractionInFam ) continue;

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

        if ( cl_4 - cl_3 > opt.lhrGapCutoff )
            ntP[ it->first ] = 1 - boost::math::cdf(X2_dist, 2*(cl_4 - cl_3) );
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

