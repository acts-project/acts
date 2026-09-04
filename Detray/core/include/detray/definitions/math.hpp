// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Algebra-plugins include(s)
#include "detray/definitions/algebra.hpp"

// SYCL include(s).
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#endif

namespace detray::math {

/// Namespace to pick up math functions from
#if defined(CL_SYCL_LANGUAGE_VERSION) || defined(SYCL_LANGUAGE_VERSION)
using namespace ::sycl;
#else
// Allow algebra plugins to define their own basic math functions
using namespace ::detray::algebra::math;
#endif  // SYCL

}  // namespace detray::math
