// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray core include(s).
#include "detray/definitions/algebra.hpp"

// Detray detector include(s)
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/telescope_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/detectors/wire_chamber_metadata.hpp"

namespace detray::test {

#define DETRAY_DEFINE_TEST_TYPES(ALGEBRA)            \
  template <detray::concepts::value T>               \
  using algebra_type = ::detray::ALGEBRA<T>;         \
                                                     \
  struct test_specialisation_name {                  \
    template <typename T>                            \
    static std::string GetName(int i) {              \
      switch (i) {                                   \
        case 0:                                      \
          return std::string(#ALGEBRA) + "<float>";  \
        case 1:                                      \
          return std::string(#ALGEBRA) + "<double>"; \
        default:                                     \
          return "unknown";                          \
      }                                              \
    }                                                \
  };

// Select algebra-plugin to compile the test with
#if DETRAY_ALGEBRA_ARRAY
static constexpr char filenames[] = "array-";
DETRAY_DEFINE_TEST_TYPES(array)

#elif DETRAY_ALGEBRA_EIGEN
static constexpr char filenames[] = "eigen-";
DETRAY_DEFINE_TEST_TYPES(eigen)

#elif DETRAY_ALGEBRA_FASTOR
static constexpr char filenames[] = "fastor-";
DETRAY_DEFINE_TEST_TYPES(fastor)

#elif DETRAY_ALGEBRA_SMATRIX
static constexpr char filenames[] = "smatrix-";
DETRAY_DEFINE_TEST_TYPES(smatrix)

#elif DETRAY_ALGEBRA_VC_AOS
static constexpr char filenames[] = "vc_aos-";
DETRAY_DEFINE_TEST_TYPES(vc_aos)

#elif DETRAY_ALGEBRA_VC_SOA
static constexpr char filenames[] = "vc_soa-";
DETRAY_DEFINE_TEST_TYPES(vc_soa)
#else
#error \
    "No algebra plugin selected for tests! Please link to one of the algebra plugins."
#endif

// Test algebra types
using algebra = algebra_type<DETRAY_CUSTOM_SCALARTYPE>;
using index = dindex_type<algebra>;
using scalar = dscalar<algebra>;
using point2 = dpoint2D<algebra>;
using point3 = dpoint3D<algebra>;
using vector2 = dvector2D<algebra>;
using vector3 = dvector3D<algebra>;
using transform3 = dtransform3D<algebra>;
template <std::size_t ROWS, std::size_t COLS>
using matrix = dmatrix<algebra, ROWS, COLS>;

// Test detector types
using default_metadata = detray::default_metadata<algebra>;
using toy_metadata = detray::toy_metadata<algebra>;
using default_telescope_metadata =
    detray::telescope_metadata<algebra, rectangle2D>;
using wire_chamber_metadata = detray::wire_chamber_metadata<algebra>;

}  // namespace detray::test
