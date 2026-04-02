#ifndef FUZZY_SCORING_H
#define FUZZY_SCORING_H

#include <algorithm>
#include <cctype>

/**
 * Nucleotide similarity (Singleton approach)
 * Returns 1.0 for an exact match, 0.0 for any mismatch.
 * This removes the intermediate fuzzy transitions (0.6/0.2).
 */
inline double nucleotideSimilarity(char a, char b)
{
    a = static_cast<char>(std::toupper(static_cast<unsigned char>(a)));
    b = static_cast<char>(std::toupper(static_cast<unsigned char>(b)));

    // Singleton logic: Membership is 1.0 ONLY at the point of identity
    if (a == b) return 1.0;

    // All other cases (Transitions/Transversions) are 0.0
    return 0.0;
}

/**
 * Fuzzy membership (Singleton-style)
 * Since similarity is now 1 or 0, match/mismatch will also be binary.
 */
inline void fuzzyMembership(double similarity,
                            double &match,
                            double &mismatch)
{
    similarity = std::clamp(similarity, 0.0, 1.0);
    
    // In singleton logic, these become indicators:
    // If similarity is 1.0, match is 1.0 and mismatch is 0.0
    // If similarity is 0.0, match is 0.0 and mismatch is 1.0
    match    = similarity;
    mismatch = 1.0 - similarity;
}

/**
 * Calculates the substitution score based on singleton membership.
 */
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