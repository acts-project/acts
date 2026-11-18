// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
  /// Material collection identifier enum supplied by the detector.
  using material_id = typename detector_t::material::id;
  /// Scalar type used by the detector's material definitions.
  using scalar_type = dscalar<typename detector_t::algebra_type>;

  static_assert(concepts::detector<detector_t>);

  /// Hash function for the homogeneous surface material (slabs and rods)
  ///
  /// Uses the material (@c std::hash) and the thickness (radius), which are
  /// compared by the equality operator of the material slab (rod).
  struct material_hash {
    /// @tparam mat_t the type of the homogeneous material (slab or rod)
    /// @param mat the material slab or rod
    /// @returns the hash value of @p mat
    template <typename mat_t>
    DETRAY_HOST std::size_t operator()(const mat_t &mat) const noexcept {
      std::size_t seed{std::hash<material<scalar_type>>{}(mat.get_material())};
      if constexpr (requires { mat.thickness(); }) {
        seed = hash_combine(seed, mat.thickness());
      } else {
        seed = hash_combine(seed, mat.radius());
      }
      return seed;
    }
  };

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
  /// @param sf_factory Factory that supplies surfaces and/or material data.
  /// @param ctx Geometry context forwarded to the decorated volume builder.
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

  /// Add the volume and the material to the detector @p det
  /// @param det Detector that receives material data.
  /// @param opt Options controlling volume construction.
  /// @param ctx Geometry context used by the decorated builder.
  /// @returns Pointer to the constructed volume.
  DETRAY_HOST
  auto build(detector_t &det, const volume_builder_options &opt,
             typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    DETRAY_VERBOSE_HOST("Build homogeneous material...");
    DETRAY_DEBUG_HOST("-> n_surfaces=" << this->surfaces().size());

    if (opt.deduplicate()) {
      DETRAY_VERBOSE_HOST("-> Deduplicate homogeneous surface material");
    }

    if constexpr (concepts::has_material_slabs<detector_t>) {
      add_material<material_id::e_material_slab>(det, opt);
    }
    if constexpr (concepts::has_material_rods<detector_t>) {
      add_material<material_id::e_material_rod>(det, opt);
    }

    m_materials.clear_all();

    DETRAY_VERBOSE_HOST(
        "Successfully built homogeneous material for volume: " << this->name());

    // Call the underlying volume builder(s) and give the volume to the
    // next decorator
    return volume_decorator<detector_t>::build(det, opt, ctx);
  }

 private:
  /// Add material of type @tparam mat_id to a detector
  ///
  /// @param det the detector
  /// @param cache the hash maps of the detector data (no deduplication if
  ///              null)
  template <material_id mat_id>
  DETRAY_HOST void add_material(detector_t &det,
                                const volume_builder_options &opt) {
    auto &local_coll = m_materials.template get<mat_id>();
    if (local_coll.empty()) {
      return;
    }

    // non-const access
    auto &det_coll = det._materials.template get<mat_id>();
    [[maybe_unused]] const auto offset{static_cast<dindex>(det_coll.size())};

    // Global index for every volume local material entry
    const std::vector<dindex> global_idx{
        this->add_objects(det_coll, local_coll, material_hash{}, opt)};

    // Update all surface links
    for (auto &sf : this->surfaces()) {
      if (sf.material().id() == mat_id) {
        sf.material().set_index(global_idx.at(sf.material().index()));
        DETRAY_DEBUG_HOST("-> material link now: " << sf.material());
      }
    }

    DETRAY_VERBOSE_HOST("-> Appended " << det_coll.size() - offset << " of "
                                       << local_coll.size()
                                       << " entries of type " << mat_id
                                       << " to detector materials");
  }

  // Material container for this volume
  typename detector_t::material_container m_materials{};
};

}  // namespace detray
