// #include <iostream>
// #include <vector>
// #include <string>
// #include <algorithm>
// #include <iomanip>
// #include <fstream>
// #include <cstdlib>
// #include <cmath>
// #include <chrono>
// #include "triangular.h"

// using namespace std;

// void smithWaterman(const string& seq1, const string& seq2, double gapPenalty) {


    

//     int m = seq1.length();
//     int n = seq2.length();

//     // DP matrix
//     vector<vector<double>> dp(m + 1, vector<double>(n + 1, 0.0));

//     double maxScore = 0;
//     int maxI = 0, maxJ = 0;

//     // Fill DP matrix
//     for (int i = 1; i <= m; i++) {
//         for (int j = 1; j <= n; j++) {

//             double matchMismatch = dp[i - 1][j - 1] +
//                 fuzzySubstitutionScore(seq1[i - 1], seq2[j - 1]);

//             double deleteGap = dp[i - 1][j] + gapPenalty;
//             double insertGap = dp[i][j - 1] + gapPenalty;

//             dp[i][j] = max({0.0, matchMismatch, deleteGap, insertGap});

//             if (dp[i][j] >= maxScore) {
//                 maxScore = dp[i][j];
//                 maxI = i;
//                 maxJ = j;
//             }
//         }
//     }

//     // ---- SAVE MATRIX TO CSV ----
//     ofstream file("matrix.csv");

//     cout<< "max:" << maxI <<", " << maxJ; 
//     // First row → column headers
//     file << "-";
//     for (int j = 0; j < n; j++) {
//         file << "," << seq2[j];
//     }
//     file << "\n";

//     // Matrix rows
//     for (int i = 0; i <= m; i++) {

//         if (i == 0)
//             file << "-";
//         else
//             file << seq1[i - 1];

//         for (int j = 0; j <= n; j++) {
//             file << "," << fixed << setprecision(2) << dp[i][j];
//         }

//         file << "\n";
//     }

//     file.close();

//     // ---- TRACEBACK ----
//     string align1 = "";
//     string align2 = "";

//     // Track the path cells
//     vector<pair<int,int>> path;

//     int i = maxI, j = maxJ;

//     double eps = 1e-6;

//     while (dp[i][j] > eps) {

//         path.push_back({i, j});

//         // Compute all three predecessor scores with bounds guards
//         double diagScore = (i > 0 && j > 0)
//             ? dp[i-1][j-1] + fuzzySubstitutionScore(seq1[i-1], seq2[j-1])
//             : -1e9;
//         double upScore   = (i > 0) ? dp[i-1][j] + gapPenalty : -1e9;
//         double leftScore = (j > 0) ? dp[i][j-1] + gapPenalty : -1e9;


//         bool fromDiag = fabs(dp[i][j] - diagScore) < eps;
//         bool fromUp   = fabs(dp[i][j] - upScore)   < eps;
//         bool fromLeft = fabs(dp[i][j] - leftScore)  < eps;

//         if (fromDiag) {
//             // Diagonal: match or mismatch
//             align1 = seq1[i-1] + align1;
//             align2 = seq2[j-1] + align2;
//             i--;
//             j--;
//         }
//         else if (fromUp) {
//             // Up: gap in seq2
//             align1 = seq1[i-1] + align1;
//             align2 = "-" + align2;
//             i--;
//         }
//         else if (fromLeft) {
//             // Left: gap in seq1
//             align1 = "-" + align1;
//             align2 = seq2[j-1] + align2;
//             j--;
//         }
//         else {
//             // Floating-point fallback: snap to the closest predecessor
//             double bestPrev = max({diagScore, upScore, leftScore});
//             if (fabs(bestPrev - diagScore) < eps) {
//                 align1 = seq1[i-1] + align1;
//                 align2 = seq2[j-1] + align2;
//                 i--; j--;
//             } else if (fabs(bestPrev - upScore) < eps) {
//                 align1 = seq1[i-1] + align1;
//                 align2 = "-" + align2;
//                 i--;
//             } else {
//                 align1 = "-" + align1;
//                 align2 = seq2[j-1] + align2;
//                 j--;
//             }
//         }
//     }


//     // ---- SAVE TRACEBACK PATH TO CSV ----
//     ofstream pathFile("traceback_path.csv");
//     pathFile << "row,col\n";
//     for (auto& p : path) {
//         pathFile << p.first << "," << p.second << "\n";
//     }
//     pathFile.close();

//     // ---- SAVE ALIGNMENT TO TXT ----
//     ofstream alignFile("alignment.txt");
//     alignFile << align1 << "\n" << align2 << "\n";
//     alignFile << fixed << setprecision(2) << maxScore << "\n";
//     alignFile.close();

//     // Output alignment in terminal
//     cout << "\nOptimal Local Alignment:\n";
//     cout << align1 << endl;
//     cout << align2 << endl;
//     cout << "Alignment Score: " << maxScore << endl;

//     // ---- AUTOMATICALLY CALL PYTHON VISUALIZATION ----
//     cout << "\nOpening DP Matrix in separate window...\n";
//     system("python visualize_matrix.py");
// }

// int main() {

//     string seq1, seq2;
//     int gapPenalty = -1;

//     cout << "Enter first DNA sequence: ";
//     cin >> seq1;

//     cout << "Enter second DNA sequence: ";
//     cin >> seq2;

//     // auto start = chrono::high_resolution_clock::now();

//     // for(int i=0; i<1000;i++)
//     // {
//     smithWaterman(seq1, seq2, gapPenalty);
//     // }

//     // auto end = chrono::high_resolution_clock::now();

//     // auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);

//     // cout << "Execution time: " << duration.count()/1000 << " microseconds\n";
   

//     return 0;
// }

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
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
// Smith-Waterman Class
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

        // DP computation
        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; j++) {

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

        // TRACEBACK
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
public:

    static void saveMatrix(const vector<vector<double>>& dp,
                           const string& seq1,
                           const string& seq2)
    {
        ofstream file("matrix.csv");

        file << "-";
        for (char c : seq2) file << "," << c;
        file << "\n";

        for (int i = 0; i <= seq1.length(); i++) {
            if (i == 0) file << "-";
            else file << seq1[i-1];

            for (int j = 0; j <= seq2.length(); j++) {
                file << "," << fixed << setprecision(2) << dp[i][j];
            }
            file << "\n";
        }
        file.close();
    }

    static void savePath(const vector<pair<int,int>>& path) {
        ofstream file("traceback_path.csv");
        file << "row,col\n";
        for (auto& p : path)
            file << p.first << "," << p.second << "\n";
        file.close();
    }

    static void saveAlignment(const AlignmentResult& res) {
        ofstream file("alignment.txt");
        file << res.align1 << "\n" << res.align2 << "\n";
        file << fixed << setprecision(2) << res.score << "\n";
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

    if (forward.score >= reverse.score) {
        best = forward;
        strand = "+";
    } else {
        best = reverse;
        strand = "-";
    }

    // Save outputs
    OutputManager::saveMatrix(best.dp, reference, read);
    OutputManager::savePath(best.path);
    OutputManager::saveAlignment(best);

    // Console output
    cout << "\nBest Alignment (Strand " << strand << "):\n";
    cout << best.align1 << endl;
    cout << best.align2 << endl;
    cout << "Score: " << best.score << endl;

    // Visualization
    system("python visualize_matrix.py");

    return 0;
}