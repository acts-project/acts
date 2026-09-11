/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/geometry/host_detector.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <string_view>

namespace traccc::io {

/// Construct a "default" / ODD detector geometry from a set of input files
///
/// It is mostly just a wrapper around Detray's more generic
/// @c detray::io::read_detector code.
///
/// @param detector The detector object to be set up
/// @param mr The memory resource to be used by the detector object
/// @param geometry_file The file containing the geometry description
/// @param material_file The file containing the material description
/// @param grid_file The file containing the detector grid description
///
void read_detector(host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view& geometry_file,
                   const std::string_view& material_file = "",
                   const std::string_view& grid_file = "");

}  // namespace traccc::io
