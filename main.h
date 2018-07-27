
#ifndef MAIN_H_
#define MAIN_H_

#define PROGRAM "lhmut"
#define VERSION "v0.1"
#define AUTHORS "Rui YE"
#define CONTACT "yerui@connect.hku.hk"
#define REMARKS "(likelihood-based mutation detector)"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <utility>
#include <getopt.h>
#include <algorithm>

using namespace std;

typedef unsigned long ulong;

class Option {
    public:
        Option():
            infileName(""),
            outfileName(""),
            baseQuaCutoff(10),
            minSupOnEachStrand(3),
            minFractionInFam(0.001),
            freqPrecision(0.0001),
            lhrGapCutoff(10.0),
            phredOffset(33),
            debug(false) {}

        ~Option() {}

        string infileName;
        string outfileName;
        int    baseQuaCutoff;
        ulong  minSupOnEachStrand;
        double minFractionInFam;
        double freqPrecision;
        double lhrGapCutoff;
        int    phredOffset;
        bool   debug;
};

typedef map<char,double> mCharDouble;
typedef map<char, ulong> mCharUlong;
typedef map<string,ulong> mStrUlong;
typedef map<string, vector<string> > mStrStrV;

typedef pair<char,ulong> pCharUlong;
typedef pair<string,ulong> pStrUlong;
typedef pair<double, set<char> > pDoubleCharSet;

typedef vector< set<char> > vCharSet;
typedef vector<string> vString;

// descending sort pair
inline bool _cmpByFirst(const pDoubleCharSet &a, const pDoubleCharSet &b) { return a.first > b.first; }
inline bool _cmpBySecond(const pCharUlong &a, const pCharUlong &b) { return a.second > b.second; }
inline bool _cmpBySecond_StrUlong(const pStrUlong &a, const pStrUlong &b) { return a.second > b.second; }


void usage();
void parseOption(int argc, char **argv, Option &opt);
void replace (string &str, const string &from, const string &to, size_t more=0 );
mStrUlong fetchInDel(string &seq, char type);
vector< pair<string, ulong> > selectInDel( const mStrUlong &m );

inline double errorRate(const char &q, const Option &opt)
{
    return pow( 10, (double)(opt.phredOffset-(int)q)/10 );
}

inline vector<double> quaToErrorRate(const string &qs, const Option &opt)
{
    vector<double> eV;
    for ( size_t i(0); i != qs.size(); i++ )  eV.push_back( errorRate(qs[i], opt) );
    return eV;
}

inline string _itoa( size_t &i )
{
    ostringstream osm;
    osm << i;
    return osm.str();
}

inline ulong countN(const string &s)
{
    ulong n(0);
    for ( size_t i(0); i != s.size(); i++ ) {
        if ( s[i] == 'N' ) n++;
    }
    return n;
}

inline bool lowQuality(const char &q, const Option &opt)
{
    return ( int(q) - opt.phredOffset < opt.baseQuaCutoff );
}

inline string reverseComplement(const string &seq)
{
    string str("");
    for ( string::const_reverse_iterator it = seq.rbegin(); it != seq.rend(); it++) {
        switch (*it) {
            case 'A': str += "T";    break;
            case 'C': str += "G";    break;
            case 'G': str += "C";    break;
            case 'T': str += "A";    break;
            case 'N':
            default:  str += "N";    break;
        }
    }

    return str;
}

inline void convertBase(string &seq, char &ref)
{
    for ( auto &i : seq ) {
        switch(i) {
            case 'a': i = 'A'; break;
            case 'c': i = 'C'; break;
            case 'g': i = 'G'; break;
            case 't': i = 'T'; break;
            case 'n': i = 'N'; break;
            case '.':
            case ',': i = ref; break;
            default:           break;
        }
    }
}

inline mCharUlong countBaseNum(const string &s)
{
    mCharUlong m;
    for ( auto &i : s ) m[i]++;
    return m;
}
#endif // MAIN_H_
