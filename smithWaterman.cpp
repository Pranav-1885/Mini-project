#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "fuzzy_scoring.h"   // or gaussian/sigmoidal/t2

using namespace std;

// ─────────────────────────────────────────────
// FASTA Reader
// ─────────────────────────────────────────────
string readFasta(const string& filename) {
ifstream file(filename);
string line, seq = "";


while (getline(file, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') continue;
    seq += line;
}
return seq;


}

// ─────────────────────────────────────────────
// FASTQ Reader
// ─────────────────────────────────────────────
struct FastqRead {
string sequence;
string quality;
};

vector<FastqRead> readFastq(const string& filename) {
ifstream file(filename);
vector<FastqRead> reads;


string h, seq, plus, qual;

while (getline(file, h)) {
    getline(file, seq);
    getline(file, plus);
    getline(file, qual);
    reads.push_back({seq, qual});
}

return reads;


}

// ─────────────────────────────────────────────
// Quality Normalization
// ─────────────────────────────────────────────
inline double qualityWeight(char q) {
double val = (q - 33) / 40.0;
return clamp(val, 0.0, 1.0);
}

// ─────────────────────────────────────────────
// Result Structure
// ─────────────────────────────────────────────
struct AlignmentResult {
string align1;
string align2;
double score;
vector<pair<int,int>> path;
vector<vector<double>> dp;
int maxI, maxJ;
};

// ─────────────────────────────────────────────
// Smith-Waterman Class (PARALLEL VERSION)
// ─────────────────────────────────────────────
class SmithWaterman {
private:
double gapPenalty;
double eps = 1e-6;

public:
SmithWaterman(double gap) : gapPenalty(gap) {}


AlignmentResult compute(const string& seq1,
                        const string& seq2,
                        const string& quality) {

    int m = seq1.length();
    int n = seq2.length();

    vector<vector<double>> dp(m + 1, vector<double>(n + 1, 0.0));

    double maxScore = 0;
    int maxI = 0, maxJ = 0;

    // Anti-diagonal parallel computation
    for (int d = 2; d <= m + n; d++) {

        int start_i = max(1, d - n);
        int end_i   = min(m, d - 1);

        for (int i = start_i; i <= end_i; i++) {

            int j = d - i;

            double base = fuzzySubstitutionScore(seq1[i - 1], seq2[j - 1]);
            double w = qualityWeight(quality[j - 1]);

            double matchMismatch = dp[i - 1][j - 1] + (base * w);
            double deleteGap = dp[i - 1][j] + gapPenalty;
            double insertGap = dp[i][j - 1] + gapPenalty;

            dp[i][j] = max({0.0, matchMismatch, deleteGap, insertGap});

            if (dp[i][j] >= maxScore) {
                maxScore = dp[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    // ───── TRACEBACK ─────
    string align1 = "", align2 = "";
    vector<pair<int,int>> path;

    int i = maxI, j = maxJ;

    while (dp[i][j] > eps) {

        path.push_back({i, j});

        double base = fuzzySubstitutionScore(seq1[i-1], seq2[j-1]);
        double w = qualityWeight(quality[j-1]);

        double diagScore = dp[i-1][j-1] + (base * w);
        double upScore   = dp[i-1][j] + gapPenalty;
        double leftScore = dp[i][j-1] + gapPenalty;

        if (fabs(dp[i][j] - diagScore) < eps) {
            align1 = seq1[i-1] + align1;
            align2 = seq2[j-1] + align2;
            i--; j--;
        }
        else if (fabs(dp[i][j] - upScore) < eps) {
            align1 = seq1[i-1] + align1;
            align2 = "-" + align2;
            i--;
        }
        else {
            align1 = "-" + align1;
            align2 = seq2[j-1] + align2;
            j--;
        }
    }

    return {align1, align2, maxScore, path, dp, maxI, maxJ};
}


};

// ─────────────────────────────────────────────
// DNA Utility Class
// ─────────────────────────────────────────────
class DNAUtils {
public:
static char complement(char c) {
switch (toupper(c)) {
case 'A': return 'T';
case 'T': return 'A';
case 'C': return 'G';
case 'G': return 'C';
case 'U': return 'A';
default: return c;
}
}


static string reverseComplement(const string& seq) {
    string rev = "";
    for (int i = seq.length() - 1; i >= 0; i--) {
        rev += complement(seq[i]);
    }
    return rev;
}


};

// ─────────────────────────────────────────────
// File Writer Class (SAM)
// ─────────────────────────────────────────────
class OutputManager {
private:
static string generateCigar(const string& a1, const string& a2, int l, int r) {
string cigar = "";
if (l > 0) cigar += to_string(l) + "S";


    int count = 0;
    char lastOp = ' ';

    for (size_t i = 0; i < a1.length(); i++) {
        char op;
        if (a1[i] == '-') op = 'I';
        else if (a2[i] == '-') op = 'D';
        else op = 'M';

        if (op == lastOp) count++;
        else {
            if (count > 0) cigar += to_string(count) + lastOp;
            lastOp = op;
            count = 1;
        }
    }
    if (count > 0) cigar += to_string(count) + lastOp;
    if (r > 0) cigar += to_string(r) + "S";

    return cigar;
}


public:
static void saveSAM(const AlignmentResult& res,
const string& refName,
int refLength,
const string& readName,
const string& usedRead,
bool isReverse)
{
ofstream file("alignment.sam");


    file << "@HD\tVN:1.6\tSO:unsorted\n";
    file << "@SQ\tSN:" << refName << "\tLN:" << refLength << "\n";

    int first_i = res.path.empty() ? 0 : res.path.back().first;
    int first_j = res.path.empty() ? 0 : res.path.back().second;

    int pos = first_i;
    int leftClip = first_j - 1;
    int rightClip = usedRead.length() - res.maxJ;

    string cigar = generateCigar(res.align1, res.align2, leftClip, rightClip);

    file << readName << "\t"
         << (isReverse ? 16 : 0) << "\t"
         << refName << "\t"
         << pos << "\t255\t"
         << cigar << "\t*\t0\t0\t"
         << usedRead << "\t*\t"
         << "AS:i:" << (int)res.score << "\n";
}


};

// ─────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────
int main() {


string fastaFile, fastqFile;

cout << "Enter FASTA file: ";
cin >> fastaFile;

cout << "Enter FASTQ file: ";
cin >> fastqFile;

string reference = readFasta(fastaFile);
auto reads = readFastq(fastqFile);

SmithWaterman sw(-1);

for (int i = 0; i < reads.size(); i++) {

    auto& r = reads[i];

    auto forward = sw.compute(reference, r.sequence, r.quality);

    string revRead = DNAUtils::reverseComplement(r.sequence);
    string revQual = string(r.quality.rbegin(), r.quality.rend());

    auto reverse = sw.compute(reference, revRead, revQual);

    bool isReverse = (reverse.score > forward.score);
    auto best = isReverse ? reverse : forward;
    string usedRead = isReverse ? revRead : r.sequence;

    OutputManager::saveSAM(best, "REF", reference.length(),
                           "READ_" + to_string(i+1),
                           usedRead, isReverse);

    cout << "\nRead " << i+1 << ":\n";
    cout << best.align1 << endl;
    cout << best.align2 << endl;
    cout << "Score: " << best.score << endl;
}

return 0;


}