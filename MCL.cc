/*
 * Modified Composite Log Likelihood
 *
 * Based on YONG CHEN, JING HUANG, et al.
 * A conditional composite likelihood ratio test with boundary constraints.
 * Biometrika (2018), 105, 1, pp. 225–232
 *
 * Written by Ye Rui (yerui@connect.hku.hk) 20181016
 *
 */
#include <string>
#include <vector>
#include <map>
#include <armadillo>

#include "MCL.h"
using namespace arma;
using namespace std;

// composite log likelihood: l_c(theta)
// mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const string         &base,
        const vector<double> &errRateV,
        const mat            &theta        )
{
    double l_c(0.0);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i];
        switch( base[i] ) {
            case 'A': l_c += log( (1-4*e/3) * theta(0) + e/3 ); break;
            case 'C': l_c += log( (1-4*e/3) * theta(1) + e/3 ); break;
            case 'G': l_c += log( (1-4*e/3) * theta(2) + e/3 ); break;
            case 'T': l_c += log( (1-4*e/3) * theta(3) + e/3 ); break;
            case 'N':                                           break;
            default: cerr << "unknown base in " << base << endl,exit(1);
        }
    }

    return l_c;
}

// composite score function: U_c(theta)
// return a column vector
//
mat composite_score (
        const string         &base,
        const vector<double> &errRateV,
        const mat            &theta        )
{
    mat U_c(4, 1, fill::zeros);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i];

        switch( base[i] ) {
            case 'A': U_c(0) += (1-4*e/3) / ( (1-4*e/3)*theta(0) + e/3 );  break;
            case 'C': U_c(1) += (1-4*e/3) / ( (1-4*e/3)*theta(1) + e/3 );  break;
            case 'G': U_c(2) += (1-4*e/3) / ( (1-4*e/3)*theta(2) + e/3 );  break;
            case 'T': U_c(3) += (1-4*e/3) / ( (1-4*e/3)*theta(3) + e/3 );  break;
            case 'N':                                                      break;
            default: cerr << "unknown base in " << base << endl, exit(1);
        }
    }

    return U_c;
}

// sensitivity matrix: H
// quaFrac is a hash table which contains base error rate with its empirical fraction on that position.
// return a 4x4 matrix, only the diagnal has values
//
mat sensitivity_matrix (
        const map<double, double> &errRateFrac,
        const mat                 &theta        )
{
    mat H(4, 4, fill::zeros);

    for (uword i(0); i != H.n_cols; i++) {
        for (map<double, double>::const_iterator it = errRateFrac.begin(); it != errRateFrac.end(); it++) {
            const double &e = it->first;
            H(i,i) += it->second * ( pow(1-4*e/3, 2) / ((1-4*e/3)*theta(i) + e/3) );
        }
    }

    return H;
}

mat sensitivity_matrix (
        const string         &base,
        const vector<double> &errRateV,
        const mat            &theta        )
{
    mat H(4, 4, fill::zeros);

    for (size_t i(0); i != base.size(); i++ ) {
        const double &e = errRateV[i];

        switch( base[i] ) {
            case 'A': H(0,0) += pow(1-4*e/3, 2) / pow( (1-4*e/3)*theta(0) + e/3, 2 );  break;
            case 'C': H(1,1) += pow(1-4*e/3, 2) / pow( (1-4*e/3)*theta(1) + e/3, 2 );  break;
            case 'G': H(2,2) += pow(1-4*e/3, 2) / pow( (1-4*e/3)*theta(2) + e/3, 2 );  break;
            case 'T': H(3,3) += pow(1-4*e/3, 2) / pow( (1-4*e/3)*theta(3) + e/3, 2 );  break;
            case 'N':                                                      break;
            default: cerr << "unknown base in " << base << endl, exit(1);
        }
    }

    return H / base.size();
}

// variability matrix: V
// return a 4x4 matrix
//
mat variability_matrix (
        const long                &N,
        const map<double, double> &errRateFrac,
        const mat                 &theta        )
{
    mat V(4, 4, fill::ones);
    double t(0.0);

    for (map<double, double>::const_iterator it = errRateFrac.begin(); it != errRateFrac.end(); it++) {
        const double &e = it->first;
        t += it->second * ( 1-4*e/3 );
    }

    return N * t*t * V;
}

mat variability_matrix (const long &N, const mat U)
{
    return U * U.t() / N;
}

// hatH_A is the inverse of the robust variance estimator
// return 4x4 matrix hatH_A
//
mat get_hatH_A (
        const mat &hatH,
        const mat &hatV        )
{
    return hatH * pinv(hatV) * hatH; // use pinv() since hatV is a singular
}

// T(theta)
// return a colume matrix
//
mat get_T_theta (
        const long &N,
        const mat  &hatH,
        const mat  &U_c_hatTheta_c,
        const mat  &theta,
        const mat  &hatTheta_c        )
{
    return pow(N, -0.5) * pinv(hatH) * U_c_hatTheta_c - pow(N, 0.5) * (theta - hatTheta_c);
}

// phi(theta)
//
double get_phi_theta (
        const double &l_c_theta,
        const double &l_c_hatTheta_c,
        const long   &N,
        const mat    &T_theta,
        const mat    &hatH,
        const mat    &U_c_hatTheta_c        )
{
    return (l_c_theta - l_c_hatTheta_c) /
        ( as_scalar(-T_theta.t() * hatH * T_theta) +
         pow(N, -1.0) * as_scalar(U_c_hatTheta_c.t() * pinv(hatH) * U_c_hatTheta_c) );
}

// l_M(theta)
//
double get_l_M_theta (
        const double    &l_c_hatTheta_c,
        const double    &phi_theta,
        const mat       &T_theta,
        const mat       &hatH_A        )
{
    return l_c_hatTheta_c - as_scalar(T_theta.t() * hatH_A * T_theta) * phi_theta;
}

// l_M(theta) from origin
//
double get_l_M_theta_origin (
        const std::string &base,
        const std::vector<double> &errRateV,
        const arma::mat &theta,
        const arma::mat &hatTheta_c,
        const std::map<double, double> &errRateFrac )
{
    size_t N = base.size();

    double l_c_theta = composite_LogLikelihood( base, errRateV, theta );
    double l_c_hatTheta_c = composite_LogLikelihood( base, errRateV, hatTheta_c );

    arma::mat U_c_hatTheta_c = composite_score( base, errRateV, hatTheta_c );
//    arma::mat hatH = sensitivity_matrix( errRateFrac, hatTheta_c );
    arma::mat hatH = sensitivity_matrix( base, errRateV, hatTheta_c );
//    arma::mat hatV = variability_matrix( N, errRateFrac, hatTheta_c );
    arma::mat hatV = variability_matrix( N, U_c_hatTheta_c );
    arma::mat hatH_A = get_hatH_A( hatH, hatV );

    arma::mat T_theta = get_T_theta( N, hatH, U_c_hatTheta_c, theta, hatTheta_c );
    double phi_theta = get_phi_theta( l_c_theta, l_c_hatTheta_c, N, T_theta, hatH, U_c_hatTheta_c );
    double l_M_theta = get_l_M_theta( l_c_hatTheta_c, phi_theta, T_theta, hatH_A );

    return l_M_theta;
}
