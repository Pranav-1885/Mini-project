// ─────────────────────────────────────────────────────────────────────────────
// fuzzy_sw_parallel.cpp
//
// Runs Smith-Waterman with all 5 fuzzy membership functions in parallel.
// Each function runs in its own thread; results are written to separate SAM
// files (alignment_singleton.sam, alignment_triangular.sam, etc.).
//
// Compile:  g++ -std=c++17 -O2 -pthread fuzzy_sw_parallel.cpp -o fuzzy_sw
// ─────────────────────────────────────────────────────────────────────────────

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <mutex>
#include <map>

// ── Membership function headers ───────────────────────────────────────────────
#include <cctype>

namespace st_singleton {
#undef FUZZY_SCORING_H
#include "singleton.h"
}
struct SingletonMF {
    static double fuzzySubstitutionScore(char a, char b) {
        return st_singleton::fuzzySubstitutionScore(a, b);
    }
};

namespace st_triangular {
#undef FUZZY_SCORING_H
#include "triangular.h"
}
struct TriangularMF {
    static double fuzzySubstitutionScore(char a, char b) {
        return st_triangular::fuzzySubstitutionScore(a, b);
    }
};

namespace st_gaussian {
#undef FUZZY_SCORING_H
#include "gaussian.h"
}
struct GaussianMF {
    static double fuzzySubstitutionScore(char a, char b) {
        return st_gaussian::fuzzySubstitutionScore(a, b);
    }
};

namespace st_sigmoidal {
#undef FUZZY_SCORING_H
#include "sigmoidal.h"
}
struct SigmoidalMF {
    static double fuzzySubstitutionScore(char a, char b) {
        return st_sigmoidal::fuzzySubstitutionScore(a, b);
    }
};

namespace st_t2_gaussian {
#undef FUZZY_SCORING_H
#include "t2_gaussian.h"
}
struct T2GaussianMF {
    static double fuzzySubstitutionScore(char a, char b) {
        return st_t2_gaussian::fuzzySubstitutionScore(a, b);
    }
};

using namespace std;

// ─────────────────────────────────────────────
// FASTA Reader
// ─────────────────────────────────────────────
string readFasta(const string& filename) {
    ifstream file(filename);
    string line, seq;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
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
        if (!getline(file, seq))  break;
        if (!getline(file, plus)) break;
        if (!getline(file, qual)) break;
        reads.push_back({seq, qual});
    }
    return reads;
}

// ─────────────────────────────────────────────
// Quality Normalization
// ─────────────────────────────────────────────
inline double qualityWeight(char q) {
    return clamp((q - 33) / 40.0, 0.0, 1.0);
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
// Smith-Waterman — templated on membership fn
// ─────────────────────────────────────────────
// MembershipFn must provide:
//   static double fuzzySubstitutionScore(char a, char b);
template<typename MembershipFn>
class SmithWaterman {
private:
    double gapPenalty;
    static constexpr double eps = 1e-6;

public:
    explicit SmithWaterman(double gap = -1.0) : gapPenalty(gap) {}

    AlignmentResult compute(const string& seq1,
                            const string& seq2,
                            const string& quality) const {
        int m = static_cast<int>(seq1.size());
        int n = static_cast<int>(seq2.size());

        vector<vector<double>> dp(m + 1, vector<double>(n + 1, 0.0));
        double maxScore = 0.0;
        int maxI = 0, maxJ = 0;

        // ── Anti-diagonal wavefront (preserves parallelism within a diagonal) ──
        for (int d = 2; d <= m + n; ++d) {
            int start_i = max(1, d - n);
            int end_i   = min(m, d - 1);

            for (int i = start_i; i <= end_i; ++i) {
                int j = d - i;
                double base = MembershipFn::fuzzySubstitutionScore(seq1[i-1], seq2[j-1]);
                double w    = qualityWeight(quality[j-1]);

                dp[i][j] = max({0.0,
                                dp[i-1][j-1] + base * w,
                                dp[i-1][j]   + gapPenalty,
                                dp[i][j-1]   + gapPenalty});

                if (dp[i][j] >= maxScore) {
                    maxScore = dp[i][j];
                    maxI = i; maxJ = j;
                }
            }
        }

        // ── Traceback ──────────────────────────────────────────────────────────
        string align1, align2;
        vector<pair<int,int>> path;
        int i = maxI, j = maxJ;

        while (i > 0 && j > 0 && dp[i][j] > eps) {
            path.push_back({i, j});
            double base  = MembershipFn::fuzzySubstitutionScore(seq1[i-1], seq2[j-1]);
            double w     = qualityWeight(quality[j-1]);
            double diag  = dp[i-1][j-1] + base * w;
            double up    = dp[i-1][j]   + gapPenalty;

            if (fabs(dp[i][j] - diag) < eps) {
                align1 = seq1[i-1] + align1;
                align2 = seq2[j-1] + align2;
                --i; --j;
            } else if (fabs(dp[i][j] - up) < eps) {
                align1 = seq1[i-1] + align1;
                align2 = '-'       + align2;
                --i;
            } else {
                align1 = '-'       + align1;
                align2 = seq2[j-1] + align2;
                --j;
            }
        }

        return {align1, align2, maxScore, path, dp, maxI, maxJ};
    }
};

// ─────────────────────────────────────────────
// DNA Utility
// ─────────────────────────────────────────────
namespace DNAUtils {
    char complement(char c) {
        switch (toupper(c)) {
            case 'A': return 'T'; case 'T': return 'A';
            case 'C': return 'G'; case 'G': return 'C';
            case 'U': return 'A'; default:  return c;
        }
    }
    string reverseComplement(const string& seq) {
        string rev;
        rev.reserve(seq.size());
        for (int i = static_cast<int>(seq.size()) - 1; i >= 0; --i)
            rev += complement(seq[i]);
        return rev;
    }
}

// ─────────────────────────────────────────────
// SAM Output
// ─────────────────────────────────────────────
namespace OutputManager {
    string generateCigar(const string& a1, const string& a2, int l, int r) {
        string cigar;
        if (l > 0) cigar += to_string(l) + 'S';
        int count = 0; char lastOp = ' ';
        for (size_t k = 0; k < a1.size(); ++k) {
            char op = (a1[k] == '-') ? 'I' : (a2[k] == '-') ? 'D' : 'M';
            if (op == lastOp) { ++count; }
            else { if (count) cigar += to_string(count) + lastOp; lastOp = op; count = 1; }
        }
        if (count) cigar += to_string(count) + lastOp;
        if (r > 0) cigar += to_string(r) + 'S';
        return cigar;
    }

    // Each membership function writes its own SAM file (no shared file handle)
    void saveSAM(const AlignmentResult& res,
                 const string& refName,
                 int refLength,
                 const string& readName,
                 const string& usedRead,
                 bool isReverse,
                 const string& filename) {
        ofstream file(filename, ios::app);   // append so all reads go to one file per MF

        // Write header only once (check if file was empty)
        if (file.tellp() == 0) {
            file << "@HD\tVN:1.6\tSO:unsorted\n";
            file << "@SQ\tSN:" << refName << "\tLN:" << refLength << "\n";
        }

        int pos       = res.path.empty() ? 0 : res.path.back().first;
        int leftClip  = res.path.empty() ? 0 : res.path.back().second - 1;
        int rightClip = static_cast<int>(usedRead.size()) - res.maxJ;
        string cigar  = generateCigar(res.align1, res.align2, leftClip, rightClip);

        file << readName   << '\t'
             << (isReverse ? 16 : 0) << '\t'
             << refName    << '\t'
             << pos        << "\t255\t"
             << cigar      << "\t*\t0\t0\t"
             << usedRead   << "\t*\t"
             << "AS:i:"    << static_cast<int>(res.score) << '\n';
    }
}

// ─────────────────────────────────────────────
// Per-thread alignment task
// ─────────────────────────────────────────────
struct ThreadResult {
    string mfName;
    vector<pair<int, double>> readScores; // (readIndex, bestScore)
};

// Runs Smith-Waterman for ALL reads using one membership function.
// Writes results to its own SAM file.  Thread-safe (no shared mutable state).
template<typename MembershipFn>
void alignAllReads(const string& mfName,
                   const string& reference,
                   const vector<FastqRead>& reads,
                   mutex& printMtx,
                   ThreadResult& out) {
    SmithWaterman<MembershipFn> sw(-1.0);
    out.mfName = mfName;

    // Truncate / create the SAM file fresh for this run
    { ofstream f("alignment_" + mfName + ".sam"); }

    for (int i = 0; i < static_cast<int>(reads.size()); ++i) {
        const auto& r = reads[i];

        auto fwd = sw.compute(reference, r.sequence, r.quality);

        string revSeq  = DNAUtils::reverseComplement(r.sequence);
        string revQual(r.quality.rbegin(), r.quality.rend());
        auto   rev     = sw.compute(reference, revSeq, revQual);

        bool   isRev   = (rev.score > fwd.score);
        auto&  best    = isRev ? rev : fwd;
        string usedRead = isRev ? revSeq : r.sequence;

        string readName = "READ_" + to_string(i + 1);
        OutputManager::saveSAM(best, "REF",
                                static_cast<int>(reference.size()),
                                readName, usedRead, isRev,
                                "alignment_" + mfName + ".sam");

        out.readScores.push_back({i + 1, best.score});

        // Print to console — lock to avoid interleaved output
        {
            lock_guard<mutex> lg(printMtx);
            cout << "[" << mfName << "] Read " << (i + 1) << ":\n"
                 << "  " << best.align1 << '\n'
                 << "  " << best.align2 << '\n'
                 << "  Score: " << best.score << '\n';
        }
    }
}

// ─────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────
int main() {
    string fastaFile, fastqFile;
    cout << "Enter FASTA file: "; cin >> fastaFile;
    cout << "Enter FASTQ file: "; cin >> fastqFile;

    string reference = readFasta(fastaFile);
    auto   reads     = readFastq(fastqFile);

    if (reference.empty()) { cerr << "Error: empty reference\n"; return 1; }
    if (reads.empty())     { cerr << "Error: no reads\n";        return 1; }

    // ── Shared resources ──────────────────────────────────────────────────────
    mutex printMtx;

    // One result slot per membership function
    vector<ThreadResult> results(5);

    // ── Launch 5 threads — one per membership function ────────────────────────
    vector<thread> threads;
    threads.reserve(5);

    threads.emplace_back(alignAllReads<SingletonMF>,
                         "singleton",  ref(reference), ref(reads),
                         ref(printMtx), ref(results[0]));

    threads.emplace_back(alignAllReads<TriangularMF>,
                         "triangular", ref(reference), ref(reads),
                         ref(printMtx), ref(results[1]));

    threads.emplace_back(alignAllReads<GaussianMF>,
                         "gaussian",   ref(reference), ref(reads),
                         ref(printMtx), ref(results[2]));

    threads.emplace_back(alignAllReads<SigmoidalMF>,
                         "sigmoidal",  ref(reference), ref(reads),
                         ref(printMtx), ref(results[3]));

    threads.emplace_back(alignAllReads<T2GaussianMF>,
                         "t2gaussian", ref(reference), ref(reads),
                         ref(printMtx), ref(results[4]));

    // ── Wait for all threads ──────────────────────────────────────────────────
    for (auto& t : threads) t.join();

    // ── Summary table ─────────────────────────────────────────────────────────
    cout << "\n╔══════════════════════════════════════════════════════╗\n";
    cout <<   "║          Score summary across membership functions   ║\n";
    cout <<   "╠══════════════════════════════════════════════════════╣\n";
    cout <<   "║  Read  │ Singleton │ Trianglr │ Gaussian │ Sigmoid │ Type-2 ║\n";
    cout <<   "╠══════════════════════════════════════════════════════╣\n";

    for (int i = 0; i < static_cast<int>(reads.size()); ++i) {
        cout << "║  " << (i + 1);
        for (auto& r : results) {
            if (i < static_cast<int>(r.readScores.size()))
                printf("  │ %8.2f", r.readScores[i].second);
        }
        cout << "  ║\n";
    }
    cout << "╚══════════════════════════════════════════════════════╝\n";
    cout << "\nSAM files written: alignment_{singleton,triangular,gaussian,"
            "sigmoidal,t2gaussian}.sam\n";

    return 0;
}