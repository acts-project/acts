// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

// System include(s).
#include <algorithm>
#include <cmath>

namespace detray::algebra {

/// Namespace to pick up common math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
namespace math {
using namespace ::sycl;
}
#else
namespace math {
using namespace std;
}
#endif  // SYCL

}  // namespace detray::algebra
