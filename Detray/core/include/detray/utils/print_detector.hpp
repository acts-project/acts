// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"

// System include(s)
#include <iomanip>
#include <sstream>
#include <string>

namespace detray::utils {

namespace detail {

/// A functor that retrieves an acceleration struct and prints it
struct accelerator_printer {
  /// Print an acceleration structure
  ///
  /// @param accel_coll collection of acceleration structs
  /// @param idx the specific grid to be checked
  /// @param id type id of the material grid collection
  template <typename accel_coll_t, typename index_t>
  DETRAY_HOST void operator()(const accel_coll_t &accel_coll, const index_t idx,
                              std::stringstream &os) const {
    os << accel_coll[idx];
  }
};

/// A functor that retrieves material and prints it
struct material_printer {
  /// Print material
  ///
  /// @param material_coll collection of material
  /// @param idx the specific grid to be checked
  /// @param id type id of the material grid collection
  template <typename material_coll_t, typename index_t>
  DETRAY_HOST void operator()(const material_coll_t &material_coll,
                              const index_t idx, std::stringstream &os) const {
    os << material_coll[idx] << std::endl;
  }
};

}  // namespace detail

/// Print basic information about the detector @param det
template <typename detector_t>
DETRAY_HOST inline std::string print_detector(
    const detector_t &det, const typename detector_t::name_map &names = {}) {
  // Gathers navigation information across navigator update calls
  std::stringstream debug_stream{};

  debug_stream << std::left << "[>] Detector " << det.name(names) << " has "
               << det.volumes().size() << " volumes." << std::endl;

  for (const auto [i, v_desc] : detray::views::enumerate(det.volumes())) {
    tracking_volume v{det, v_desc};

    debug_stream << "[>>] Volume " << v.name(names) << std::endl;
    debug_stream << v << std::endl;

    debug_stream << "[>>>] Acceleration Structures:" << std::endl;
    const auto acc_link = v_desc.accel_link();
    for (std::size_t j = 0u; j < acc_link.size(); ++j) {
      // An acceleration data structure link was set, but is invalid
      if (!acc_link[j].is_invalid_id() && !acc_link[j].is_invalid_index()) {
        debug_stream << "Surfaces registered in '" << acc_link[j].id()
                     << "':" << std::endl;
        det.accelerator_store().template visit<detail::accelerator_printer>(
            acc_link[j], debug_stream);
      }
    }

    debug_stream << "[>>>] Surfaces:" << std::endl;
    for (const auto sf_desc : v.template surfaces<>()) {
      geometry::surface sf{det, sf_desc};
      debug_stream << sf << std::endl;

      // Check the surface material, if present
      if (sf.has_material()) {
        debug_stream << "[>>>>] Surface material: ";
        sf.template visit_material<detail::material_printer>(debug_stream);
      }
    }

    debug_stream << std::endl;
  }

  return debug_stream.str();
}

}  // namespace detray::utils
