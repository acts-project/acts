/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/definitions/qualifiers.hpp"

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

// System include(s).
#include <cmath>

namespace traccc::math {

#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
using ::sycl::abs;
using ::sycl::acos;
using ::sycl::asin;
using ::sycl::atan;
using ::sycl::atan2;
using ::sycl::cos;
using ::sycl::exp;
using ::sycl::fabs;
using ::sycl::floor;
using ::sycl::fmod;
using ::sycl::log;
using ::sycl::max;
using ::sycl::min;
using ::sycl::pow;
using ::sycl::sin;
using ::sycl::sqrt;
using ::sycl::tan;
#else
using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cos;
using std::exp;
using std::fabs;
using std::floor;
using std::fmod;
using std::log;
using std::max;
using std::min;
using std::pow;
using std::sin;
using std::sqrt;
using std::tan;
#endif

/**
 * @brief Perform IEEE-754 division, even if fast math is enabled.
 *
 * @returns x divided by y
 */
template <typename T>
TRACCC_HOST_DEVICE inline __attribute__((always_inline)) T div_ieee754(T x,
                                                                       T y) {
    static_assert(std::is_same_v<T, double> || std::is_same_v<T, float>);
#ifdef __CUDA_ARCH__
    if constexpr (std::is_same_v<T, double>) {
        return __ddiv_rn(x, y);
    } else {
        return __fdiv_rn(x, y);
    }
#else
    return x / y;
#endif
}
}  // namespace traccc::math
