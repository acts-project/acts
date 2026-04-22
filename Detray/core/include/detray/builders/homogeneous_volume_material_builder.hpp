// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <iostream>
#include <map>
#include <memory>

namespace detray {

/// @brief Build a volume with homogeneous volume material.
///
/// Decorator class to a volume builder that adds the material data to the
/// volume while building the volume.
template <typename detector_t>
class homogeneous_volume_material_builder final
    : public volume_decorator<detector_t> {
 public:
  using scalar_type = dscalar<typename detector_t::algebra_type>;

  /// @param vol_builder volume builder that should be decorated with volume
  /// material
  DETRAY_HOST
  explicit homogeneous_volume_material_builder(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : volume_decorator<detector_t>(std::move(vol_builder)) {}

  /// Add all necessary components for the volume material
  ///
  /// @param mat the material parameters
  DETRAY_HOST
  void set_material(material<scalar_type> mat) { m_volume_material = mat; }

  /// Add the volume and the material to the detector @param det
  DETRAY_HOST
  auto build(detector_t &det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    // Call the underlying volume builder(s)
    typename detector_t::volume_type *vol =
        volume_decorator<detector_t>::build(det, ctx);

    // Nothing left to do
    if (m_volume_material == detray::vacuum<scalar_type>{}) {
      DETRAY_WARN_HOST("Volume " << this->vol_index()
                                 << " has vacuum material");

      return vol;
    }

    constexpr auto material_id{detector_t::material::id::e_raw_material};

    // Update the volume material link
    dindex coll_size{det.material_store().template size<material_id>()};
    vol->set_material(material_id, coll_size);

    // Append the material
    det._materials.template push_back<material_id>(m_volume_material);

    // Give the volume to the next decorator
    return vol;
  }

 private:
  // Material for this volume
  material<scalar_type> m_volume_material{};
};

}  // namespace detray
