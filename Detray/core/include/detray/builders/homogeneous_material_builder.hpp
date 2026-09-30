// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/builders/detail/material_deduplication.hpp"
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

  static_assert(concepts::detector<detector_t>);

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

  /// Toggles whether material that is identical to material already present
  /// in the detector is shared instead of being copied
  DETRAY_HOST
  void deduplicate_material(bool toggle) override {
    m_deduplicate = toggle;
    volume_decorator<detector_t>::deduplicate_material(toggle);
  }

  /// @returns whether identical material is shared between surfaces
  DETRAY_HOST
  bool deduplicate_material() const { return m_deduplicate; }

  /// Add the volume and the material to the detector @param det
  DETRAY_HOST
  auto build(detector_t &det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    DETRAY_VERBOSE_HOST("Build homogeneous material...");

    DETRAY_DEBUG_HOST("-> n_surfaces=" << this->surfaces().size());

    if (m_deduplicate) {
      add_deduplicated_material(det);
    } else {
      add_material(det);
    }

    DETRAY_VERBOSE_HOST(
        "Successfully built homogeneous material for volume: " << this->name());

    // Call the underlying volume builder(s) and give the volume to the
    // next decorator
    return volume_decorator<detector_t>::build(det, ctx);
  }

 private:
  /// Append the material of the volume to the detector @param det and shift
  /// the surface material links accordingly
  DETRAY_HOST
  void add_material(detector_t &det) {
    const auto &material = det.material_store();

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
    det._materials.append(std::move(m_materials));
    m_materials.clear_all();
  }

  /// Add only the material of the volume to the detector @param det that is
  /// not yet present there and link the surfaces to the existing entries
  /// otherwise
  DETRAY_HOST
  void add_deduplicated_material(detector_t &det) {
    DETRAY_VERBOSE_HOST("-> Deduplicate homogeneous material");

    if constexpr (concepts::has_material_slabs<detector_t>) {
      deduplicate<material_id::e_material_slab>(det);
    }
    if constexpr (concepts::has_material_rods<detector_t>) {
      deduplicate<material_id::e_material_rod>(det);
    }

    // Add remaining material types, if any
    det._materials.append(std::move(m_materials));
    m_materials.clear_all();
  }

  /// Deduplicate the material of type @tparam mat_id against the material in
  /// the detector @param det
  template <material_id mat_id>
  DETRAY_HOST void deduplicate(detector_t &det) {
    auto &local_coll = m_materials.template get<mat_id>();
    if (local_coll.empty()) {
      return;
    }

    auto &det_coll = det._materials.template get<mat_id>();
    [[maybe_unused]] const std::size_t n_before{det_coll.size()};

    // Global index for every volume local material entry
    detail::homogeneous_material_lookup lookup{det_coll};
    std::vector<dindex> global_idx;
    global_idx.reserve(local_coll.size());
    for (const auto &mat : local_coll) {
      global_idx.push_back(lookup.insert(mat));
    }

    for (auto &sf : this->surfaces()) {
      if (sf.material().id() == mat_id) {
        sf.material().set_index(global_idx.at(sf.material().index()));
        DETRAY_DEBUG_HOST("-> sf=" << sf
                                   << ": material now: " << sf.material());
      }
    }

    DETRAY_DEBUG_HOST("-> Appended " << det_coll.size() - n_before << " of "
                                     << local_coll.size() << " entries of type "
                                     << mat_id << " to detector materials");
    local_coll.clear();
  }

  /// Whether to share identical material between surfaces
  bool m_deduplicate{false};
  // Material container for this volume
  typename detector_t::material_container m_materials{};
};

}  // namespace detray
