#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <algorithm>

#include "compass_search.h"
#include "main.h"
#include "likelihood.h"

vector<pDoubleCharSet> llh_genotype(const string &s, const string &q, const Option &opt)
{
    set<char> ntS;
    vector<pDoubleCharSet> llhV;

    // check frequent of alleles
    mCharUlong fr;
    double depth( s.size() );

    for ( size_t i(0); i != depth; i++ ) {
        if ( lowQuality(q[i], opt) || s[i] == 'N' )  continue;
        fr[ s[i] ]++;
    }

    if ( fr.empty() ) {
        return llhV;
    }

    // sort by frequent
    vector<pCharUlong> ntV( fr.begin(), fr.end() );
/*
    if ( opt.debug ) {
        for ( auto &p : ntV ) {
            cout << p.first << ':' << p.second << ' ';
        }
        cout << '\n';
    }
*/
    // only one allele
    if ( ntV.size() == 1 ) {
        if ( ntV[0].second >= opt.minSupOnEachStrand ) {
            ntS.insert( ntV[0].first );
            llhV.push_back( make_pair(1.0, ntS) );
        }

        return llhV;
    }

    if ( ntV.size() > 1 )
        sort( ntV.begin(), ntV.end(), _cmpBySecond ); // descending sort

    // delete pair<allele, supportNum> which has too few support reads or too small fraction
    while ( ntV.size() ) {
        if ( ntV.back().second < opt.minSupOnEachStrand || ntV.back().second / depth < opt.minFractionInFam )
            ntV.pop_back();
        else
            break;
    }

    // none allele remain
    if ( ntV.empty() ) {
        return llhV;
    }

    // only one allele
    if ( ntV.size() == 1 ) {
        if ( ntV[0].second >= opt.minSupOnEachStrand ) {
            ntS.insert( ntV[0].first );
            llhV.push_back( make_pair(1.0, ntS) );
        }

        return llhV;
    }

    // two or more alleles
    vector<double> errV = quaToErrorRate(q, opt);

    if ( ntV.size() > 1 ) {
        vector<pCharUlong> tmp(ntV.begin()+1, ntV.end());

        llhV.push_back( maxLogLikelihood(s,errV,ntV,opt,1) ); // 1st
        llhV.push_back( maxLogLikelihood(s,errV,tmp,opt,1) ); // 2nd

        llhV.push_back( maxLogLikelihood(s,errV,ntV,opt,2) );


        if ( ntV.size() > 2 ) {
            llhV.push_back( maxLogLikelihood(s,errV,tmp,opt,2) ); // 2nd & 3rd

            tmp.erase(tmp.begin());
            llhV.push_back( maxLogLikelihood(s,errV,tmp,opt,1) ); // 3rd

            tmp = ntV;
            tmp.erase(tmp.begin()+1);
            llhV.push_back( maxLogLikelihood(s,errV,tmp,opt,2) ); // 1st & 3rd


            llhV.push_back(  maxLogLikelihood(s,errV,ntV,opt,3) );
        }

    }

    for ( auto & p : llhV ) {
        if ( p.second.size() == 3 ) p.first -= 2 * opt.lhrGapCutoff;
        if ( p.second.size() == 2 ) p.first -=     opt.lhrGapCutoff;
    }

    if ( llhV.size() > 1 )
        sort(llhV.begin(), llhV.end(), _cmpByFirst); // desending sort

    // four alleles is really rare, so do not consider it here
    return llhV;
}

pDoubleCharSet maxLogLikelihood(const string &s, const vector<double> &e, const vector<pCharUlong> &v, const Option &opt, const int mode)
{
    double llh(0.0);
    set<char> cSet;

    if ( mode == 1 ) {
        char a( v[0].first );
        for ( size_t i(0); i != s.size(); i++ ) {
            if ( s[i] == a )
                llh += log( 1 - e[i] );
            else
                llh += log( e[i]/3 );
        }

        cSet.insert(a);
        return make_pair(llh, cSet);
    }
    else if ( mode == 2 ) {
        //initial value
        char a( v[0].first ), b( v[1].first );
        double total( v[0].second + v[1].second );
        double af( v[0].second / total);
        double bf( 1 - af );

        for ( size_t i(0); i != s.size(); i++ ) {
            if ( s[i] == 'N' ) continue;

            s[i] == a ? ( llh += log( (1 - e[i] * 4/3) * af + e[i]/3 ) )
                : s[i] == b ? ( llh += log( (1 - e[i] * 4/3) * bf + e[i]/3 ) )
                :             ( llh += log( e[i] / 3 )                       )
                ;
        }

        double mllh = llh; // max log likelihood
        double oldaf = af;
        long step(1);

        if (opt.debug) cout << "~org af bf llh: " << af << ' ' << bf << ' ' << llh << endl;

        // calculate max likelihood step by step
        while (1) {
            af += opt.freqPrecision * step;
            if ( af > 1 ) af = 1;
            if ( af < 0 ) af = 0;
            bf = 1 - af;

            llh = 0;
            for ( size_t i(0); i != s.size(); i++ ) {
                if ( s[i] == 'N' ) continue;

                s[i] == a ? ( llh += log( (1 - e[i] * 4/3) * af + e[i]/3 ) )
                    : s[i] == b ? ( llh += log( (1 - e[i] * 4/3) * bf + e[i]/3 ) )
                    :             ( llh += log( e[i]/3 )                         )
                    ;
            }

            if (opt.debug) cout << "~loop af bf step llh: " << af << ' ' << bf << ' ' << step << ' ' << llh << endl;

            if ( llh > mllh ) {
                mllh = llh;
                step *= 2;  // the direction is right, increase step length
                oldaf = af;
                if ( af == 1 ) step = -1;
                if ( af == 0 ) step = 1;
            }
            else {
                af = oldaf;
                if ( step == 1 )
                    step = -1;
                else if ( step == -1) {
                    bf = 1 - af;
                    break;
                }
                else
                    step = (step > 0 ? 1 : -1);
            }
        }

        if ( af != 0 )   cSet.insert(a);
        if ( bf != 0 )   cSet.insert(b);

        return make_pair(mllh, cSet);
    }
    else if ( mode == 3 ) {
        //
        double total(0.0);
        //        mCharDouble mBaseFreq;
        for ( size_t i(0); i != 3; i++ ) {
            total += v[i].second;
        }
        //
        // although 3 alleles, only two free variables, maximize the function of two variables.

        double delta(0.01);
        double delta_tol(0.00001);
        double fx;
        int k;
        int k_max(20000);
        int m = 2;

        double *x;
        double *x0;
        double *lower;
        double *upper;

        lower = new double[m];
        upper = new double[m];
        lower[0] = 0.0;
        upper[0] = 1.0;
        lower[1] = 0.0;
        upper[1] = 1.0;

        x0 = new double[m];
        if (opt.debug) cout << "  Test COMPASS_SEARCH with the function.\n";

        x0[0] = v[0].second/total;
        x0[1] = v[1].second/total;

        if (opt.debug) {
            r8vec_print ( m, x0, "  Initial point X0:" );
            cout << "\n";
            cout << "  F(X0) = " << minus_llh_3nt( m, x0, v, s, e, lower, upper, 1.0 ) << "\n";
        }

        x = compass_search ( minus_llh_3nt, m, x0, v, s, e, lower, upper, 1.0, delta_tol, delta, k_max, fx, k );

        if (opt.debug) {
            r8vec_print ( m, x, "  Estimated minimizer X1:" );
            cout << "\n";
            cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
        }

        // first two allele
        for ( int i(0); i != 2; i++ ) {
            if ( x[i] >= opt.freqPrecision  )   cSet.insert(v[i].first);
        }

        // the third allele
        if ( 1 - x[0] - x[1] >= opt.freqPrecision )   cSet.insert(v[2].first);

        delete [] lower;
        delete [] upper;
        delete [] x;
        delete [] x0;

        return make_pair(-fx,cSet); // -fx is the biggest llh
    }
    else {
        cerr << "only support mode 1,2,3" << endl, exit(1);
    }
}

// minimizing -llh equals to maximizing llh
double minus_llh_3nt( int m, double x[],
        const vector<pCharUlong> &v,
        const string &s, const vector<double> &e,
        double lower[], double upper[],
        double sumBound )
{
    double llh(0.0);

    for ( int i(0); i != m; i++ ) {
        if ( x[i] < lower[i] ) x[i] = lower[i];
        if ( x[i] > upper[i] ) x[i] = upper[i];
    }

    if ( x[0] + x[1] > sumBound ) x[1] = sumBound - x[0];

    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) continue;

        s[i] == v[0].first ? ( llh += log( (1 - e[i] * 4/3) * x[0] + e[i]/3 ) )
            : s[i] == v[1].first ? ( llh += log( (1 - e[i] * 4/3) * x[1] + e[i]/3 ) )
            : s[i] == v[2].first ? ( llh += log( (1 - e[i] * 4/3) * ( 1-x[0]-x[1] ) + e[i]/3 ) )
            :                      ( llh += log( e[i] / 3 )                       )
            ;
    }

    return -llh;
}
