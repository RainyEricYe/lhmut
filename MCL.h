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

#ifndef _MCL_H_
#define _MCL_H_

// composite log likelihood: l_c(theta)
// arma::mat theta is a column vector which has 4 elements for A, C, G, T, respectively.
//
double composite_LogLikelihood (
        const std::string    &base,
        const std::vector<double> &errRateV,
        const arma::mat      &theta );

// composite score function: U_c(theta)
// return a column vector
//
arma::mat composite_score (
        const std::string    &base,
        const std::vector<double> &errRateV,
        const arma::mat      &theta );

// sensitivity matrix: H
// quaFrac is a hash table which contains base error rate with its empirical fraction on that position.
// return a 4x4 matrix, only the diagnal has values
//
arma::mat sensitivity_matrix (
        const std::map<double, double> &errRateFrac,
        const arma::mat                &theta );

arma::mat sensitivity_matrix (
        const std::string         &base,
        const std::vector<double> &errRateV,
        const arma::mat           &theta        );

// variability matrix: V
// return a 4x4 matrix
//
arma::mat variability_matrix (
        const long                     &N,
        const std::map<double, double> &errRateFrac,
        const arma::mat                &theta );

// variability matrix: V
// return a 4x4 matrix based on composite_score U
//
arma::mat variability_matrix (
        const long       &N,
        const arma::mat  U  );

// hatH_A is the inverse of the robust variance estimator
// return 4x4 matrix hatH_A
//
arma::mat get_hatH_A (
        const arma::mat &hatH,
        const arma::mat &hatV );

// T(theta)
// return a colume matrix
//
arma::mat get_T_theta (
        const long      &N,
        const arma::mat &hatH,
        const arma::mat &U_c_hatTheta_c,
        const arma::mat &theta,
        const arma::mat &hatTheta_c );

// phi(theta)
//
double get_phi_theta (
        const double    &l_c_theta,
        const double    &l_c_hatTheta_c,
        const long      &N,
        const arma::mat &T_theta,
        const arma::mat &hatH,
        const arma::mat &U_c_hatTheta_c );

// l_M(theta)
//
double get_l_M_theta (
        const double    &l_c_hatTheta_c,
        const double    &phi_theta,
        const arma::mat &T_theta,
        const arma::mat &hatH_A );

// l_M(theta) from origin
//
double get_l_M_theta_origin (
        const std::string &base,
        const std::vector<double> &errRateV,
        const arma::mat &theta,
        const arma::mat &hatTheta_c,
        const std::map<double, double> &errRateFrac );

#endif // _MCL_H_

