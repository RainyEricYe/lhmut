#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

struct fn_data {
    std::string         base;
    std::vector<double> errRateV;
};

mCharDouble llh_genotype(const string &s, const string &q, const Option &opt);

#endif // LIKELIHOOD_H_
