// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

/// This file contains the definitions of the linear algebra types to be used in
/// the tutorials, which in this case will be the std::array based plugin.
/// The includes for the different algebra plugins can be found in
/// @c detray/plugins/algebra/

// Detray core include(s).
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"

// Detray detector include(s)
#include "detray/detectors/default_metadata.hpp"
#include "detray/detectors/telescope_metadata.hpp"
#include "detray/detectors/toy_metadata.hpp"

namespace detray::tutorial {

using algebra_t = detray::array<DETRAY_CUSTOM_SCALARTYPE>;
using scalar = dscalar<algebra_t>;
using point2 = detray::dpoint2D<algebra_t>;
using point3 = detray::dpoint3D<algebra_t>;
using vector3 = detray::dvector3D<algebra_t>;
using transform3 = detray::dtransform3D<algebra_t>;

// Test detector types
using default_metadata = detray::default_metadata<algebra_t>;
using toy_metadata = detray::toy_metadata<algebra_t>;
using default_telescope_metadata =
    detray::telescope_metadata<algebra_t, rectangle2D>;

}  // namespace detray::tutorial
