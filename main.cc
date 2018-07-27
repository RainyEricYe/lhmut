
#include "main.h"

void usage() {
    cout << "Program: " << PROGRAM << "  " << REMARKS << "\n"
        "Version: " << VERSION << "\n"
        "Authors: " << AUTHORS << "\n"
        "Contact: " << CONTACT << "\n\n"

        "Options: " << PROGRAM << " -i in.pileup -o outfile\n\n"
        "    -i [s]     input pileup file\n"
        "    -o [s]     output mutation file\n"

        "    -q [i]     base quality cutoff [10]\n"
        "    -s [i]     min support num on each strand [3]\n"
        "    -f [f]     min fraction of alterative allele in a read family [0.001]\n"

        "    -e [f]     precision of allele frequency while calculate likelihood ratio [0.0001]\n"
        "    -g [f]     gap between likelihood ratios of major and rest genotypes [10.0]\n"
        "    -x [i]     Encoding offset for phred quality scores [33]\n"

        "    -d         debug mode [false]\n"
        "    -h         help\n"
        "    -v         version\n"
        "\n";
}

void parseOption(int argc, char **argv, Option &opt) {

    if ( argc < 2 ) usage(), exit(1);

    int c;
    while ( (c=getopt(argc,argv,"i:o:q:s:f:e:g:x:dvh")) != -1 ) {
        switch (c) {
            case 'i': opt.infileName    = optarg;          break;
            case 'o': opt.outfileName   = optarg;          break;

            case 'q': opt.baseQuaCutoff      = atoi(optarg);    break;
            case 's': opt.minSupOnEachStrand = atoi(optarg);    break;
            case 'f': opt.minFractionInFam   = atof(optarg);    break;

            case 'e': opt.freqPrecision = atof(optarg);         break;
            case 'g': opt.lhrGapCutoff  = atof(optarg);         break;
            case 'x': opt.phredOffset   = atoi(optarg);         break;

            case 'd': opt.debug = true;                             break;
            case 'v': cerr << VERSION << endl;                  exit(1);
            case 'h':
            default:  usage();                                  exit(1);
        }
    }

    return;
}

void replace (string &str, const string &from, const string &to, size_t more )
{
    size_t pos;
    size_t offset = 0;
    const size_t increment = to.size();

    while ((pos = str.find(from, offset)) != string::npos) {
        str.replace(pos, from.size()+more, to);
        offset = pos + increment;
    }
}

mStrUlong fetchInDel(string &s, char type)
{
    mStrUlong m;
//    cout << "old seq: " << s << endl;
    size_t p(0);

    while ( (p=s.find(type,p)) != string::npos ) {
        string len("");
        size_t offset;
        size_t length;

        for ( size_t i(p+1); i != s.size(); i++ ) {
            if ( isdigit( s[i] ) ) {
                len += s[i];
            }
            else {
                offset = i;
                break;
            }
        }

        length = atoi(len.c_str());
        string indel( s.begin()+offset, s.begin()+offset+length );
        m[indel]++;

        s.replace(p, offset-p+length, "");
  //      cout << offset << ' '<< p << ' ' << length << ' ' << offset-p+length << ' ' << indel
//           << "\nnew seq: " << s << endl;

    }

//    for ( auto &i:m ) cout << i.first << ' ' << i.second << ' ';
  //  cout << endl;

    return m;
}

vector<pStrUlong> selectInDel( const mStrUlong &m )
{
    vector<pStrUlong> v( m.begin(), m.end() );

    if ( v.size() > 1 )
        sort( v.begin(), v.end(), _cmpBySecond_StrUlong );

    for ( vector<pStrUlong>::iterator it = v.begin(); it != v.end(); it++ ) {
        if ( countN(it->first) / (double)it->first.size() > 0.5 ) v.erase(it);
    }

    return v;
}
