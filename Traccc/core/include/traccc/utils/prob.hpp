/** TRACCC library, part of the ACTS project (R&D line)
 *
 * This file includes code from the ROOT (https://github.com/root-project/root)
 * and Cephes Library (http://www.netlib.org/cephes)
 *
 * ROOT is licensed under the GNU Lesser General Public License v2.1
 *
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

/* Logarithm of gamma function */
/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */

template <typename scalar_t>
struct log_gamma {

    TRACCC_HOST_DEVICE static scalar_t A(unsigned int i) {
        switch (i) {
            case 0:
                return A0;
            case 1:
                return A1;
            case 2:
                return A2;
            case 3:
                return A3;
            case 4:
                return A4;
        }
        assert(false);
        return 0.f;
    }

    TRACCC_HOST_DEVICE static scalar_t B(unsigned int i) {
        switch (i) {
            case 0:
                return B0;
            case 1:
                return B1;
            case 2:
                return B2;
            case 3:
                return B3;
            case 4:
                return B4;
            case 5:
                return B5;
        }
        assert(false);
        return 0.f;
    }

    TRACCC_HOST_DEVICE static scalar_t C(unsigned int i) {
        switch (i) {
            case 0:
                return C0;
            case 1:
                return C1;
            case 2:
                return C2;
            case 3:
                return C3;
            case 4:
                return C4;
            case 5:
                return C5;
        }
        assert(false);
        return 0.f;
    }

    static constexpr scalar_t kMAXLGM = static_cast<scalar_t>(2.556348e305);
    static constexpr scalar_t kMACHEP = 1.11022302462515654042363166809e-16f;
    static constexpr scalar_t kMAXLOG = 709.782712893383973096206318587f;
    static constexpr scalar_t kBig = 4.503599627370496e15f;
    static constexpr scalar_t kBiginv = 2.22044604925031308085e-16f;
    /* log( sqrt( 2*pi ) ) */
    static constexpr scalar_t LS2PI = 0.91893853320467274178f;

    private:
    static constexpr scalar_t A0 = 8.11614167470508450300E-4f;
    static constexpr scalar_t A1 = -5.95061904284301438324E-4f;
    static constexpr scalar_t A2 = 7.93650340457716943945E-4f;
    static constexpr scalar_t A3 = -2.77777777730099687205E-3f;
    static constexpr scalar_t A4 = 8.33333333333331927722E-2f;

    static constexpr scalar_t B0 = -1.37825152569120859100E3f;
    static constexpr scalar_t B1 = -3.88016315134637840924E4f;
    static constexpr scalar_t B2 = -3.31612992738871184744E5f;
    static constexpr scalar_t B3 = -1.16237097492762307383E6f;
    static constexpr scalar_t B4 = -1.72173700820839662146E6f;
    static constexpr scalar_t B5 = -8.53555664245765465627E5f;

    static constexpr scalar_t C0 = -3.51815701436523470549E2f;
    static constexpr scalar_t C1 = -1.70642106651881159223E4f;
    static constexpr scalar_t C2 = -2.20528590553854454839E5f;
    static constexpr scalar_t C3 = -1.13933444367982507207E6f;
    static constexpr scalar_t C4 = -2.53252307177582951285E6f;
    static constexpr scalar_t C5 = -2.01889141433532773231E6f;
};

// Forward declarations of functions
template <typename scalar_t>
TRACCC_HOST_DEVICE scalar_t chisquared_cdf_c(const scalar_t x,
                                             const scalar_t r);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igam(const scalar_t a, const scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igam_impl(const scalar_t a,
                                             const scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igamc(const scalar_t a, const scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igamc_impl(const scalar_t a,
                                              const scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t lgam(scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t lgam_impl(scalar_t x);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_A(const scalar_t x,
                                                    const unsigned int N);
template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_B(const scalar_t x,
                                                    const unsigned int N);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_C(const scalar_t x,
                                                    const unsigned int N);

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t prob(const scalar_t chi2,
                                        const scalar_t ndf) {
    if (ndf <= 0)
        return 0;  // Set CL to zero in case ndf<=0

    if (chi2 <= 0) {
        if (chi2 < 0)
            return 0;
        else
            return 1;
    }

    const scalar_t ret_val = chisquared_cdf_c(chi2, ndf);
    assert(ret_val >= 0.f && ret_val <= 1.f);
    return ret_val;
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t chisquared_cdf_c(const scalar_t x,
                                                    const scalar_t r) {
    scalar_t retval = igamc(0.5f * r, 0.5f * x);
    return static_cast<scalar_t>(retval);
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igam(const scalar_t a, const scalar_t x) {

    // LM: for negative values returns 1.0 instead of zero
    // This is correct if a is a negative integer since Gamma(-n) = +/- inf
    if (a <= 0)
        return 1.0;

    if (x <= 0)
        return 0.0;

    // for (x > 1) && (x > a where a > 0)
    if ((x > 1.0) && (x > a))
        return (1.0 - igamc_impl<scalar_t>(a, x));

    // for (0 < x < 1) || ( x < a where a > 0)
    return igam_impl<scalar_t>(a, x);
}

// for (0 < x < 1) || ( x < a where a > 0)
template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igam_impl(const scalar_t a,
                                             const scalar_t x) {
    scalar_t ans, ax, c, r;

    /* Compute  x**a * exp(-x) / gamma(a)  */
    ax = a * math::log(x) - x - lgam(a);
    if (ax < -log_gamma<scalar_t>::kMAXLOG)
        return (0.0f);

    ax = std::exp(ax);

    /* power series */
    r = a;
    c = 1.0f;
    ans = 1.0f;

    do {
        r += 1.0f;
        c *= x / r;
        ans += c;
    } while (c / ans > log_gamma<scalar_t>::kMACHEP);

    return (ans * ax / a);
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igamc(const scalar_t a, const scalar_t x) {

    // LM: for negative values returns 0.0
    // This is correct if a is a negative integer since Gamma(-n) = +/- inf
    if (a <= 0.f)
        return 0.0f;

    if (x <= 0.f)
        return 1.0f;

    // for (0 < x < 1) || (x < a where a > 0)
    if ((x < 1.0f) || (x < a))
        return (1.0f - igam_impl<scalar_t>(a, x));

    // for (x > 1) && (x > a where a > 0)
    return igamc_impl<scalar_t>(a, x);
}

// for (x > 1) && (x > a where a > 0)
template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t igamc_impl(const scalar_t a,
                                              const scalar_t x) {
    scalar_t ans, ax, c, yc, r, t, y, z;
    scalar_t pk, pkm1, pkm2, qk, qkm1, qkm2;

    ax = a * std::log(x) - x - lgam(a);
    if (ax < -log_gamma<scalar_t>::kMAXLOG)
        return (0.0f);

    ax = std::exp(ax);

    /* continued fraction */
    y = 1.0f - a;
    z = x + y + 1.0f;
    c = 0.0f;
    pkm2 = 1.0f;
    qkm2 = x;
    pkm1 = x + 1.0f;
    qkm1 = z * x;
    ans = pkm1 / qkm1;

    do {
        c += 1.0f;
        y += 1.0f;
        z += 2.0f;
        yc = y * c;
        pk = pkm1 * z - pkm2 * yc;
        qk = qkm1 * z - qkm2 * yc;
        if (qk != 0.f) {
            r = pk / qk;
            t = std::abs((ans - r) / r);
            ans = r;
        } else {
            t = 1.0f;
        }
        pkm2 = pkm1;
        pkm1 = pk;
        qkm2 = qkm1;
        qkm1 = qk;
        if (std::abs(pk) > log_gamma<scalar_t>::kBig) {
            pkm2 *= log_gamma<scalar_t>::kBiginv;
            pkm1 *= log_gamma<scalar_t>::kBiginv;
            qkm2 *= log_gamma<scalar_t>::kBiginv;
            qkm1 *= log_gamma<scalar_t>::kBiginv;
        }
    } while (t > log_gamma<scalar_t>::kMACHEP);

    return (ans * ax);
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t lgam(scalar_t x) {
    scalar_t p, q, u, w, z;
    int i;

    [[maybe_unused]] int sgngam = 1;

    if (x >= std::numeric_limits<scalar_t>::infinity())
        return (std::numeric_limits<scalar_t>::infinity());

    if (x < -34.0f) {
        q = -x;

        // For x > 34
        w = lgam_impl<scalar_t>(q);
        p = std::floor(q);
        if (p == q)  //_unur_FP_same(p,q)
            return (std::numeric_limits<scalar_t>::infinity());
        i = static_cast<int>(p);
        if ((i & 1) == 0)
            sgngam = -1;
        else
            sgngam = 1;
        z = q - p;
        if (z > 0.5f) {
            p += 1.0f;
            z = p - q;
        }
        z = q * std::sin(constant<scalar_t>::pi * z);
        if (z == 0)
            return (std::numeric_limits<scalar_t>::infinity());
        /* z = log(ROOT::Math::Pi()) - log( z ) - w;*/
        z = std::log(constant<scalar_t>::pi) - math::log(z) - w;
        return (z);
    }

    if (x < 13.0f) {
        z = 1.0f;
        p = 0.0f;
        u = x;
        while (u >= 3.0f) {
            p -= 1.0f;
            u = x + p;
            z *= u;
        }
        while (u < 2.0f) {
            if (u == 0.f)
                return (std::numeric_limits<scalar_t>::infinity());
            z /= u;
            p += 1.0f;
            u = x + p;
        }
        if (z < 0.0f) {
            sgngam = -1;
            z = -z;
        } else
            sgngam = 1;
        if (u == static_cast<scalar_t>(2.0))
            return (std::log(z));
        p -= 2.0f;
        x = x + p;
        p = x * Polynomialeval_B(x, 5) / Polynomialeval_C(x, 6);
        return (std::log(z) + p);
    }

    return lgam_impl<scalar_t>(x);
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t lgam_impl(scalar_t x) {
    scalar_t p, q;
    scalar_t sgngam = 1.f;

    if (x > log_gamma<scalar_t>::kMAXLGM)
        return (sgngam * std::numeric_limits<scalar_t>::infinity());

    q = (x - 0.5f) * std::log(x) - x + log_gamma<scalar_t>::LS2PI;
    if (x > 1.0e8f)
        return (q);

    p = 1.0f / (x * x);
    if (x >= 1000.0f)
        q +=
            ((7.9365079365079365079365e-4f * p - 2.7777777777777777777778e-3f) *
                 p +
             0.0833333333333333333333f) /
            x;
    else
        q += Polynomialeval_A(p, 4) / x;
    return (q);
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_A(const scalar_t x,
                                                    const unsigned int N) {

    if (N == 0)
        return log_gamma<scalar_t>::A(0);
    else {
        scalar_t pom = log_gamma<scalar_t>::A(0);
        for (unsigned int i = 1; i <= N; i++)
            pom = pom * x + log_gamma<scalar_t>::A(i);
        return pom;
    }
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_B(const scalar_t x,
                                                    const unsigned int N) {

    if (N == 0)
        return log_gamma<scalar_t>::B(0);
    else {
        scalar_t pom = log_gamma<scalar_t>::B(0);
        for (unsigned int i = 1; i <= N; i++)
            pom = pom * x + log_gamma<scalar_t>::B(i);
        return pom;
    }
}

template <typename scalar_t>
TRACCC_HOST_DEVICE inline scalar_t Polynomialeval_C(const scalar_t x,
                                                    const unsigned int N) {

    if (N == 0)
        return log_gamma<scalar_t>::C(0);
    else {
        scalar_t pom = x + log_gamma<scalar_t>::C(0);
        for (unsigned int i = 1; i < N; i++)
            pom = pom * x + log_gamma<scalar_t>::C(i);
        return pom;
    }
}

}  // namespace traccc
