#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "triangular.h"

using namespace std;

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

    AlignmentResult compute(const string& seq1, const string& seq2) {

        int m = seq1.length();
        int n = seq2.length();

        vector<vector<double>> dp(m + 1, vector<double>(n + 1, 0.0));

        double maxScore = 0;
        int maxI = 0, maxJ = 0;

        //  Anti-diagonal parallel computation
        for (int d = 2; d <= m + n; d++) {

            int start_i = max(1, d - n);
            int end_i   = min(m, d - 1);

            for (int i = start_i; i <= end_i; i++) {

                int j = d - i;

                double matchMismatch = dp[i - 1][j - 1] +
                    fuzzySubstitutionScore(seq1[i - 1], seq2[j - 1]);

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

            double diagScore = (i > 0 && j > 0)
                ? dp[i-1][j-1] + fuzzySubstitutionScore(seq1[i-1], seq2[j-1])
                : -1e9;

            double upScore   = (i > 0) ? dp[i-1][j] + gapPenalty : -1e9;
            double leftScore = (j > 0) ? dp[i][j-1] + gapPenalty : -1e9;

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
// File Writer Class
// ─────────────────────────────────────────────
class OutputManager {
private:
    // Helper to generate a CIGAR string from aligned strings and soft clips
    static string generateCigar(const string& align1, const string& align2, 
                                int leftClip, int rightClip) {
        string cigar = "";
        if (leftClip > 0) cigar += to_string(leftClip) + "S";
        
        int count = 0;
        char lastOp = ' ';
        
        for (size_t i = 0; i < align1.length(); ++i) {
            char op;
            if (align1[i] == '-') op = 'I';       // Insertion to reference
            else if (align2[i] == '-') op = 'D'; // Deletion from reference
            else op = 'M';                       // Match or Mismatch
            
            if (op == lastOp) {
                count++;
            } else {
                if (count > 0) cigar += to_string(count) + lastOp;
                lastOp = op;
                count = 1;
            }
        }
        if (count > 0) cigar += to_string(count) + lastOp;
        if (rightClip > 0) cigar += to_string(rightClip) + "S";
        
        return cigar;
    }

public:
    static void saveSAM(const AlignmentResult& res,
                        const string& refName,
                        int refLength,
                        const string& readName,
                        const string& usedRead,
                        const string& originalRead,
                        bool isReverse)
    {
        ofstream file("alignment.sam");

        // 1. Write Header
        file << "@HD\tVN:1.6\tSO:unsorted\n";
        file << "@SQ\tSN:" << refName << "\tLN:" << refLength << "\n";

        // Determine POS and soft-clips from the traceback path
        // path stores coordinates. The last element is the start of alignment (since traceback goes backwards)
        // first_i is 1-based index on reference, first_j is 1-based index on read.
        int first_i = res.path.empty() ? 0 : res.path.back().first;
        int first_j = res.path.empty() ? 0 : res.path.back().second;

        // POS (1-based leftmost mapping position)
        int pos = first_i; 

        // Calculate soft clips
        int leftClip = first_j - 1;
        int rightClip = usedRead.length() - res.maxJ;

        // Generate CIGAR string
        string cigar = generateCigar(res.align1, res.align2, leftClip, rightClip);

        // Required SAM Fields
        string qname = readName;
        int flag = isReverse ? 16 : 0;
        string rname = refName;
        int mapq = 255;          // 255 = mapping quality not available
        string rnext = "*";      // No mate/next segment
        int pnext = 0;           // No mate
        int tlen = 0;            // No template length
        
        // SEQ is the sequence mapped to the forward strand of the reference.
        // Since usedRead was aligned against the forward reference, it is already oriented properly.
        string seq = usedRead;
        string qual = "*";       // Quality unavailable

        // Optional tags: AS is the alignment score.
        string optTag = "AS:i:" + to_string(static_cast<int>(res.score));

        // Write Alignment Line
        file << qname << "\t"
             << flag << "\t"
             << rname << "\t"
             << pos << "\t"
             << mapq << "\t"
             << cigar << "\t"
             << rnext << "\t"
             << pnext << "\t"
             << tlen << "\t"
             << seq << "\t"
             << qual << "\t"
             << optTag << "\n";

        file.close();
    }
};

// ─────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────
int main() {

    string reference, read;
    double gapPenalty = -1;

    cout << "Enter reference sequence: ";
    cin >> reference;

    cout << "Enter read sequence: ";
    cin >> read;

    SmithWaterman sw(gapPenalty);

    // Forward
    AlignmentResult forward = sw.compute(reference, read);

    // Reverse
    string revRead = DNAUtils::reverseComplement(read);
    AlignmentResult reverse = sw.compute(reference, revRead);

    // Choose best
    AlignmentResult best;
    string strand;
    string usedRead;

    if (forward.score >= reverse.score) {
        best = forward;
        strand = "+";
        usedRead = read;
    } else {
        best = reverse;
        strand = "-";
        usedRead = revRead;
    }

    // Save outputs
    bool isReverse = (strand == "-");
    OutputManager::saveSAM(best, "REF", reference.length(), "READ1", usedRead, read, isReverse);

    // Console output
    cout << "\nBest Alignment (Strand " << strand << "):\n";
    cout << best.align1 << endl;
    cout << best.align2 << endl;
    cout << "Score: " << best.score << endl;
    cout << "SAM alignment exported to alignment.sam" << endl;

    return 0;
}