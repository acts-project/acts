/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include <detray/definitions/algebra.hpp>

// System include(s)
#include <cstdint>

namespace traccc {

using measurement_id_type = unsigned int;
using particle_id = std::uint64_t;
using geometry_id = std::uint64_t;
using channel_id = unsigned int;

// Default algebra type
using default_algebra = detray::array<DETRAY_CUSTOM_SCALARTYPE>;

using detray::algebra::array::operator*;
using detray::algebra::array::operator-;
using detray::algebra::array::operator+;

using scalar = detray::dscalar<default_algebra>;
using point2 = detray::dpoint2D<default_algebra>;
using vector2 = point2;
using variance2 = point2;
using point3 = detray::dpoint3D<default_algebra>;
using vector3 = detray::dvector3D<default_algebra>;
using variance3 = point3;
using transform3 = detray::dtransform3D<default_algebra>;

namespace getter = detray::getter;
namespace vector = detray::vector;
namespace matrix = detray::matrix;

// Pull in the print operator definitions for the algebra types
using detray::algebra::operator<<;

}  // namespace traccc
