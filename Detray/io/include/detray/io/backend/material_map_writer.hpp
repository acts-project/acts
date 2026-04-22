// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/backend/detail/grid_writer.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/backend/homogeneous_material_writer.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/material/material_slab.hpp"

// System include(s)
#include <string_view>

namespace detray::io {

/// @brief Material maps writer backend
///
/// Fills a material @c detector_grids_payload from a @c detector instance
class material_map_writer : public detail::grid_writer {
  using base_type = detail::grid_writer;
  using grid_writer_t = base_type;
  using mat_writer_t = homogeneous_material_writer;

 public:
  /// Tag the writer as "material_map"
  static constexpr std::string_view tag = "material_maps";

  /// Payload type that the reader processes
  using payload_type =
      detector_grids_payload<surface_material_payload, io::material_id>;

  /// Same constructors for this class as for base_type
  using base_type::base_type;

  /// Convert the header information into its payload
  template <class detector_t>
  static auto header_to_payload(const detector_t& det,
                                const std::string_view det_name) {
    return grid_writer_t::header_to_payload(tag, det.material_store(),
                                            det_name);
  }

  /// Convert the material description of a detector @param det into its io
  /// payload
  template <class detector_t>
  static payload_type to_payload(const detector_t& det,
                                 const typename detector_t::name_map&) {
    using algebra_t = typename detector_t::algebra_type;
    using material_t = material_slab<typename detector_t::scalar_type>;

    payload_type grids_data;

    for (const auto& vol_desc : det.volumes()) {
      // Volume local surface indices
      dindex offset{dindex_invalid};

      /// Check if a surface has a metrial map
      auto vol = tracking_volume{det, vol_desc};
      for (const auto& sf_desc : vol.surfaces()) {
        if (sf_desc.index() < offset) {
          offset = sf_desc.index();
        }

        const auto& mat_link = sf_desc.material();
        // Don't look at empty links
        if (mat_link.is_invalid() ||
            mat_link.id() == detector_t::material::id::e_none) {
          continue;
        }

        // Do not attempt to convert homogeneous material
        if constexpr (detray::concepts::has_material_slabs<detector_t>) {
          if (mat_link.id() == detector_t::material::id::e_material_slab) {
            continue;
          }
        }
        if constexpr (detray::concepts::has_material_rods<detector_t>) {
          if (mat_link.id() == detector_t::material::id::e_material_rod) {
            continue;
          }
        }

        // How to convert a material slab in the grid
        auto mat_converter = [&sf_desc](const material_t& mat) {
          return mat_writer_t::to_payload<algebra_t>(mat, sf_desc.index());
        };

        // Generate the payload
        grid_writer_t::to_payload(det.material_store(), mat_link,
                                  vol_desc.index(), sf_desc.index() - offset,
                                  grids_data, mat_converter);
      }
    }

    return grids_data;
  }
};

}  // namespace detray::io
