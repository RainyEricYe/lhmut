/*
 */

#include "main.h"
#include "likelihood.h"

int main( int argc, char **argv )
{
    Option opt;
    parseOption(argc, argv, opt);

    ifstream inf(opt.infileName.c_str());
    ofstream outf(opt.outfileName.c_str());

    if ( !inf.is_open() )  cerr << "open error: " << opt.infileName << endl, exit(1);
    if ( !outf.is_open() ) cerr << "open error: " << opt.outfileName << endl, exit(1);

    outf << "Chrom\tTemplate\tPos\tDepths\tMuts\t"
        "Tcount\tCcount\tGcount\tAcount\t"
        "inscount\tdelcount\tNcount" << endl;

    string line;
    while ( getline(inf, line) ) {
        string chr, baseStr, quaStr;
        ulong pos(0), depth(0);
        char ref;

        istringstream itm(line);
        itm >> chr >> pos >> ref >> depth >> baseStr >> quaStr;

        replace(baseStr, "^", "", 1); // delete lable after ^, be careful of ^$
        replace(baseStr, "$", "");
        convertBase(baseStr, ref); // to upper & , . to ref
        mCharUlong mBaseNum = countBaseNum(baseStr);



//        mStrUlong mInsNum, mDelNum;
        if ( baseStr.find("+") != string::npos ) {
//            mInsNum = fetchInDel(baseStr, '+');
            continue;
        }

        if ( baseStr.find("-") != string::npos ) {
  //          mDelNum = fetchInDel(baseStr, '-');
            continue;
        }

        depth = baseStr.size();
        if ( depth != quaStr.size() ) {
            cerr << "parse error: length of bases and qualities are different "
                << chr << ' ' << pos << endl;
            continue;
        }

        vector<pDoubleCharSet> llhV = llh_genotype(baseStr, quaStr, opt);

//        if ( llhV[0].second.size() <= 1 ) {
//            continue;
  //      }

        cout << "~~" << chr << '\t' << pos << '\t' << ref << '\t' << depth << '\t';
        string bs = "TCGAN";
        for ( auto &b : bs ) {
            mCharUlong::const_iterator i = mBaseNum.find(b);
            if ( i != mBaseNum.end() )
                cout << i->second << '\t';
            else
                cout << "0\t";
        }

        for ( auto &p : llhV ) {
            cout << p.first << ':';
            for ( auto &c : p.second ) cout << c << ' ';
        }

        if ( llhV.size() > 1 && opt.debug ) {
            vector<pCharUlong> pv(mBaseNum.begin(), mBaseNum.end());
            sort(pv.begin(), pv.end(), _cmpBySecond);

            cout << '\t' << pv.begin()->first << '\t' << 1 - pv.begin()->second/(double)depth;
            cout << '\t' << llhV[0].first - llhV[1].first;
        }

        cout << endl;
    }

    inf.close();
    outf.close();

    return 0;
}
