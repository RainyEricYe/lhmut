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

    outf <<
        "Chrom\tTemplate\tPos\tDepths\tMuts\t"
        "Tcount\tCcount\tGcount\tAcount\t"
        "inscount\tdelcount\tNcount"
        << endl;

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

        mStrUlong mInsNum = fetchInDel(baseStr, '+');
        mStrUlong mDelNum = fetchInDel(baseStr, '-');

        mCharUlong mBaseNum = countBaseNum(baseStr);

        depth = baseStr.size();
        if ( depth != quaStr.size() ) {
            cerr << "parse error: length of bases and qualities are different "
                << chr << ' ' << pos << endl;
            continue;
        }

        vector<pDoubleCharSet> llhV = llh_genotype(baseStr, quaStr, opt);


        ulong mutN(0);
        for ( auto &p : mBaseNum ) {
            if ( p.first == ref || p.first == '*' ) continue;  // not ref, ins, del
            mutN += p.second;
        }

        cout << chr << '\t' << pos << '\t' << ref << '\t'
             << depth << '\t' << mutN << '\t';

        string bs = "TCGA";
        for ( auto &b : bs ) {
            mCharUlong::const_iterator i = mBaseNum.find(b);
            cout << ( i == mBaseNum.end() ? 0 : i->second ) << '\t';
        }

        vector<pStrUlong> vInsNum = selectInDel(mInsNum);
        vector<pStrUlong> vDelNum = selectInDel(mDelNum);

        cout << ( vInsNum.empty() ? 0 : vInsNum[0].second ) << '\t';
        cout << ( vDelNum.empty() ? 0 : vDelNum[0].second ) << '\t';

        mCharUlong::const_iterator nit = mBaseNum.find('N');
        cout << ( nit == mBaseNum.end() ? 0 : nit->second );


        // detail info below

        for ( auto &p : vInsNum ) cout << "\tI+" << p.first << ':' << p.second;
        for ( auto &p : vDelNum ) cout << "\tD-" << p.first << ':' << p.second;

        for ( auto &p : llhV ) {
            cout << '\t' << p.first << ':';
            for ( auto &c : p.second ) cout << c;
        }

        if ( llhV.size() > 1 && opt.debug ) {
            vector<pCharUlong> pv(mBaseNum.begin(), mBaseNum.end());
            sort(pv.begin(), pv.end(), _cmpBySecond);

            cout << '\t' << pv.begin()->first << '\t' << 1 - pv.begin()->second/(double)depth;
            cout << '\t' << llhV[0].first - llhV[1].first << "\tdetail";
        }

        cout << endl;
    }

    inf.close();
    outf.close();

    return 0;
}
