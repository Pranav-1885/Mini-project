#ifndef FUZZY_SCORING_H
#define FUZZY_SCORING_H

#include <algorithm>
#include <cctype>
#include <cmath>

// Nucleotide similarity (crisp → [0,1])
inline double nucleotideSimilarity(char a, char b)
{
    a = std::toupper(a);
    b = std::toupper(b);

    if (a == b) return 1.0;

    // Transition pairs (purine↔purine or pyrimidine↔pyrimidine — biologically closer)
    if ((a == 'A' && b == 'G') || (a == 'G' && b == 'A') ||
        (a == 'C' && b == 'T') || (a == 'T' && b == 'C') ||
        (a == 'C' && b == 'U') || (a == 'U' && b == 'C') ||
        (a == 'T' && b == 'U') || (a == 'U' && b == 'T'))
        return 0.6;

    // Transversion (purine↔pyrimidine — less similar)
    return 0.2;
}

// Fuzzy membership (gaussian-style, dual output)
inline void fuzzyMembership(double similarity,
                            double &match,
                            double &mismatch)
{
    similarity = std::clamp(similarity, 0.0, 1.0);
    const double c = 1.0;      // center (perfect match)
    const double sigmaU = 0.35; // upper Gaussian spread
    const double sigmaL = 0.25; // lower Gaussian spread

    // Upper Gaussian membership
    double upper = exp(-pow(similarity - c, 2) / (2 * sigmaU * sigmaU));

    // Lower Gaussian membership
    double lower = exp(-pow(similarity - c, 2) / (2 * sigmaL * sigmaL));

    // Type-2 membership (average of upper and lower)
    match = (upper + lower) / 2.0;

    mismatch = 1.0 - match;
}

inline double fuzzySubstitutionScore(char a, char b)
{
    constexpr double MATCH_REWARD     =  2.0;
    constexpr double MISMATCH_PENALTY =  1.0;  

    double similarity = nucleotideSimilarity(a, b);

    double match, mismatch;
    fuzzyMembership(similarity, match, mismatch);

    return (MATCH_REWARD * match) - (MISMATCH_PENALTY * mismatch);
}

#endif