#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

vector<pDoubleCharSet> llh_genotype(const string &, const string &, const Option &);

pDoubleCharSet maxLogLikelihood(const string &, const vector<double> &, const vector<pCharUlong> &, const Option &, const int);

double minus_llh_3nt( int m, double x[],
        const vector<pCharUlong> &v,
        const string &s, const vector<double> &e,
        double lower[], double upper[],
        double sumBound );


#endif // LIKELIHOOD_H_
