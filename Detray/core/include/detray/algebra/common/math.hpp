// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

// System include(s).
#include <algorithm>
#include <cmath>

// NOLINTBEGIN
namespace detray::algebra {

/// Namespace to pick up common math functions from
///
/// The using-directive makes all functions of the backend (std or sycl)
/// available. A plugin can add overloads for its SIMD types to this namespace.
/// A qualified call (e.g. @c algebra::math::sqrt) then finds a name that is
/// declared in this namespace first, and it does not look through the
/// using-directive any more. The using-declarations below keep the backend
/// functions in the overload set of the functions that the plugins overload.
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math {
using namespace ::sycl;

using ::sycl::abs;
using ::sycl::acos;
using ::sycl::asin;
using ::sycl::atan;
using ::sycl::atan2;
using ::sycl::atanh;
using ::sycl::ceil;
using ::sycl::copysign;
using ::sycl::cos;
using ::sycl::cosh;
using ::sycl::exp;
using ::sycl::fabs;
using ::sycl::floor;
using ::sycl::fma;
using ::sycl::hypot;
using ::sycl::log;
using ::sycl::log10;
using ::sycl::max;
using ::sycl::min;
using ::sycl::pow;
using ::sycl::signbit;
using ::sycl::sin;
using ::sycl::sinh;
using ::sycl::sqrt;
using ::sycl::tan;
using ::sycl::tanh;
}  // namespace math
#else
namespace math {
using namespace std;

using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::atanh;
using std::ceil;
using std::copysign;
using std::cos;
using std::cosh;
using std::exp;
using std::fabs;
using std::floor;
using std::fma;
using std::hypot;
using std::log;
using std::log10;
using std::max;
using std::min;
using std::pow;
using std::signbit;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;
}  // namespace math
#endif  // SYCL

}  // namespace detray::algebra
// NOLINTEND
