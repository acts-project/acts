// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray core include(s).
#include "algebra/utils/data_generator.hpp"
#include "detray/definitions/algebra.hpp"

// Detray detector include(s)
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/detectors/wire_chamber_metadata.hpp"

namespace detray::benchmarks {

#define DETRAY_DEFINE_BENCHMARK_NAME(ALGEBRA) \
  static std::string plugin_name(#ALGEBRA);

// Select algebra-plugin to compile the test with
#if DETRAY_ALGEBRA_ARRAY
using algebra = detray::array<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(array)

#elif DETRAY_ALGEBRA_EIGEN
using algebra = detray::eigen<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(eigen)

#elif DETRAY_ALGEBRA_FASTOR
using algebra = detray::fastor<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(fastor)

#elif DETRAY_ALGEBRA_SMATRIX
using algebra = detray::smatrix<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(smatrix)

#elif DETRAY_ALGEBRA_VC_AOS || DETRAY_ALGEBRA_VC_SOA

#if DETRAY_ALGEBRA_VC_AOS
using algebra = detray::vc_aos<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(vc_aos)
#endif

#if DETRAY_ALGEBRA_VC_SOA
using algebra = detray::vc_soa<DETRAY_CUSTOM_SCALARTYPE>;
DETRAY_DEFINE_BENCHMARK_NAME(vc_soa)
#endif

#else
#error \
    "No algebra plugin selected for benchmarks! Please link to one of the algebra plugins."
#endif

// Test algebra types
using scalar = dscalar<algebra>;
using point2 = dpoint2D<algebra>;
using point3 = dpoint3D<algebra>;
using vector3 = dvector3D<algebra>;
using transform3 = dtransform3D<algebra>;
template <std::size_t ROWS, std::size_t COLS>
using matrix = dmatrix<algebra, ROWS, COLS>;

// Test detector types
using default_metadata = detray::default_metadata<algebra>;
using toy_metadata = detray::toy_metadata<algebra>;
using wire_chamber_metadata = detray::wire_chamber_metadata<algebra>;

}  // namespace detray::benchmarks
