/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/utils/prob.hpp"

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <Math/ProbFuncMathCore.h>
#endif  // TRACCC_HAVE_ROOT

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

// Test gamma function for pvalue evaulation
TEST(pvalue, compare_with_root) {

    // zero chisquare -> pvalue = 1
    EXPECT_FLOAT_EQ(prob(0.f, 0.1f), 1.f);
    EXPECT_FLOAT_EQ(prob(0.f, 1.f), 1.f);
    EXPECT_FLOAT_EQ(prob(0.f, 10.f), 1.f);

    // Negative chi square -> pvalue = 0
    EXPECT_FLOAT_EQ(prob(-1.f, 12.f), 0.f);
    EXPECT_FLOAT_EQ(prob(-194.f, 12.f), 0.f);

    // Non-positive NDF -> pvalue = 0
    EXPECT_FLOAT_EQ(prob(3.f, -17.4f), 0.f);
    EXPECT_FLOAT_EQ(prob(27.f, 0.f), 0.f);

    // Compare with the outputs of ROOT::Math::Prob(chi2, ndf)
    EXPECT_FLOAT_EQ(prob(3.f, 5.f), 0.69998584f);
    EXPECT_FLOAT_EQ(prob(17.4f, 6.26f), 0.0094291369f);

#ifdef TRACCC_HAVE_ROOT
    EXPECT_FLOAT_EQ(prob(3.f, 5.f),
                    static_cast<float>(ROOT::Math::chisquared_cdf_c(3., 5.)));
    EXPECT_FLOAT_EQ(
        prob(17.4f, 6.26f),
        static_cast<float>(ROOT::Math::chisquared_cdf_c(17.4, 6.26)));

    for (int i_n = 0; i_n < 500; i_n++) {
        const traccc::scalar ndf = static_cast<traccc::scalar>(i_n) * 0.1f;
        for (int i_c = 0; i_c < 100; i_c++) {

            const traccc::scalar chi2 = static_cast<traccc::scalar>(i_c) * 0.1f;
            const traccc::scalar pval = prob(ndf, chi2);
            const traccc::scalar pval_root =
                static_cast<traccc::scalar>(ROOT::Math::chisquared_cdf_c(
                    static_cast<double>(ndf), static_cast<double>(chi2)));

            ASSERT_NEAR(pval - pval_root, 0.f, pval_root * 1e-4f);
            ASSERT_GE(pval, 0.0f);
            ASSERT_LE(pval, 1.0f);
        }
    }
#endif  // TRACCC_HAVE_ROOT
}
