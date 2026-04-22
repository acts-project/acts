// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/homogeneous_material_generator.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/core/concepts.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <vector>

namespace detray {

/// @brief Build a volume containing surfaces with material.
///
/// Decorator class to a volume builder that adds the material data to the
/// surfaces while building the volume.
template <typename detector_t>
class homogeneous_material_builder final : public volume_decorator<detector_t> {
 public:
  using material_id = typename detector_t::material::id;
  using scalar_type = dscalar<typename detector_t::algebra_type>;

  /// @param vol_builder volume builder that should be decorated with material
  DETRAY_HOST
  explicit homogeneous_material_builder(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : volume_decorator<detector_t>(std::move(vol_builder)) {
    DETRAY_VERBOSE_HOST(
        "Add hom. material builder to volume: " << this->name());

    static_assert(concepts::has_material_slabs<detector_t> ||
                      concepts::has_material_rods<detector_t>,
                  "No homogeneous surface material in detector type");
  }

  /// Overwrite, to add material in addition to surfaces (only if surfaces are
  /// present in the factory, otherwise only add material)
  /// @{
  DETRAY_HOST
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) override {
    DETRAY_VERBOSE_HOST("Add [material] surface factory:");

    // If the factory also holds surface data, call base volume builder
    volume_decorator<detector_t>::add_surfaces(sf_factory, ctx);

    // Add material
    auto mat_factory =
        std::dynamic_pointer_cast<homogeneous_material_factory<detector_t>>(
            sf_factory);
    if (mat_factory) {
      DETRAY_VERBOSE_HOST("-> found decoration: " << DETRAY_TYPENAME(
                              homogeneous_material_factory<detector_t>));
      (*mat_factory)(this->surfaces(), m_materials);
      return;
    }
    auto mat_generator =
        std::dynamic_pointer_cast<homogeneous_material_generator<detector_t>>(
            sf_factory);
    if (mat_generator) {
      DETRAY_VERBOSE_HOST("-> found decoration: " << DETRAY_TYPENAME(
                              homogeneous_material_generator<detector_t>));
      (*mat_generator)(this->surfaces(), m_materials);
      return;
    }

    if (!mat_factory && !mat_generator) {
      DETRAY_VERBOSE_HOST("No material found in this surface factory");
      DETRAY_VERBOSE_HOST("-> Built non-material surfaces");
    }
  }
  /// @}

  /// Add the volume and the material to the detector @param det
  DETRAY_HOST
  auto build(detector_t &det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    DETRAY_VERBOSE_HOST("Build homogeneous material...");

    const auto &material = det.material_store();

    DETRAY_DEBUG_HOST("-> n_surfaces=" << this->surfaces().size());

    // Update the surface material links and shift them according to the
    // number of material slabs/rods that were in the detector previously
    for (auto &sf : this->surfaces()) {
      DETRAY_DEBUG_HOST("-> sf=" << sf);
      DETRAY_DEBUG_HOST("  -> material_id=" << sf.material().id());
      if constexpr (concepts::has_material_slabs<detector_t>) {
        if (sf.material().id() == material_id::e_material_slab) {
          dindex offset =
              material.template size<material_id::e_material_slab>();
          DETRAY_DEBUG_HOST("-> update material slab offset: " << offset);
          sf.update_material(offset);
          DETRAY_DEBUG_HOST("-> material now: " << sf.material());
        }

        DETRAY_DEBUG_HOST(
            "-> Appending "
            << m_materials.template size<material_id::e_material_slab>()
            << " slabs into detector materials");
      }
      if constexpr (concepts::has_material_rods<detector_t>) {
        if (sf.material().id() == material_id::e_material_rod) {
          DETRAY_DEBUG_HOST(
              "-> update material rod offset: "
              << material.template size<material_id::e_material_rod>());
          sf.update_material(
              material.template size<material_id::e_material_rod>());
        }

        DETRAY_DEBUG_HOST(
            "-> Appending "
            << m_materials.template size<material_id::e_material_rod>()
            << " rods into detector materials");
      }
    }

    // Add material to the detector

    if constexpr (types::contains<typename detector_t::material,
                                  material_rod<scalar_type>>) {
      DETRAY_DEBUG_HOST(
          "-> Appending "
          << m_materials.template size<material_id::e_material_rod>()
          << " rods into detector materials");
    }
    det._materials.append(std::move(m_materials));
    m_materials.clear_all();

    DETRAY_VERBOSE_HOST(
        "Successfully built homogeneous material for volume: " << this->name());

    // Call the underlying volume builder(s) and give the volume to the
    // next decorator
    return volume_decorator<detector_t>::build(det, ctx);
  }

 private:
  // Material container for this volume
  typename detector_t::material_container m_materials{};
};

}  // namespace detray
