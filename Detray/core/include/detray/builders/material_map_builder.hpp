// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/builders/bin_fillers.hpp"
#include "detray/builders/detail/material_deduplication.hpp"
#include "detray/builders/material_map_factory.hpp"
#include "detray/builders/material_map_generator.hpp"
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/core/concepts.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/material/material_map.hpp"

// System include(s)
#include <cassert>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace detray {

namespace detail {

template <typename material_t>
struct add_sf_material_map;

}  // namespace detail

/// @brief Build the material maps for a given volume.
///
/// Decorator class to a volume builder that adds material maps to either
/// surfaces or volumes
template <concepts::detector detector_t, std::size_t DIM = 2u,
          typename mat_map_factory_t =
              material_grid_factory<typename detector_t::algebra_type>>
class material_map_builder final : public volume_decorator<detector_t> {
  using material_t = typename detector_t::material;

 public:
  using scalar_type = dscalar<typename detector_t::algebra_type>;
  using detector_type = detector_t;
  using value_type = material_slab<scalar_type>;
  using bin_data_type =
      typename fill_by_bin::template bin_data<DIM, value_type>;

  /// @param vol_builder volume builder that should be decorated with material
  /// maps
  DETRAY_HOST
  explicit material_map_builder(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : volume_decorator<detector_t>(std::move(vol_builder)) {
    DETRAY_VERBOSE_HOST("Add material map builder to volume: " << this->name());
  }

  /// @returns the raw materials and their local bin indices that are
  /// currently staged in the builder
  DETRAY_HOST
  auto data() const -> const std::map<dindex, std::vector<bin_data_type>>& {
    return m_bin_data;
  }

  /// Overwrite, to add material maps in addition to surfaces
  /// @{
  DETRAY_HOST
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) override {
    DETRAY_VERBOSE_HOST("Add [material] surface factory:");

    // If the factory also holds surface data, call base volume builder
    volume_decorator<detector_t>::add_surfaces(sf_factory, ctx);

    // Add material and bin data
    using sf_factory_t = material_map_factory<detector_t, axis::multi_bin<DIM>>;
    auto mat_factory = std::dynamic_pointer_cast<sf_factory_t>(sf_factory);
    if (mat_factory) {
      DETRAY_VERBOSE_HOST(
          "-> found decoration: " << DETRAY_TYPENAME(sf_factory_t));
      (*mat_factory)(this->surfaces(), m_bin_data, m_n_bins, m_axis_spans);
    }
    auto mat_generator =
        std::dynamic_pointer_cast<material_map_generator<detector_t>>(
            sf_factory);
    if (mat_generator) {
      DETRAY_VERBOSE_HOST("-> found decoration: " << DETRAY_TYPENAME(
                              material_map_generator<detector_t>));
      (*mat_generator)(this->surfaces(), this->masks(), m_bin_data, m_n_bins);
      return;
    }

    if (!mat_factory && !mat_generator) {
      DETRAY_VERBOSE_HOST("No material found in this surface factory");
      DETRAY_VERBOSE_HOST("-> Built non-material surfaces");
    }
  }
  /// @}

  /// Toggles whether material maps that are identical to a map already
  /// present in the detector are shared instead of being copied
  DETRAY_HOST
  void deduplicate_material(bool toggle) override {
    m_deduplicate = toggle;
    volume_decorator<detector_t>::deduplicate_material(toggle);
  }

  /// @returns whether identical material maps are shared between surfaces
  DETRAY_HOST
  bool deduplicate_material() const { return m_deduplicate; }

  /// Add the volume and the material maps to the detector @param det
  DETRAY_HOST
  auto build(detector_t& det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type* override {
    DETRAY_VERBOSE_HOST("Build material maps...");

    // Ensure the material links are correct BEFORE the surfaces are built
    // and potentially added to an acceleration data structure
    add_material_maps(det);

    DETRAY_DEBUG_HOST(
        "-> Let underlying builders construct the volume using correct "
        "material links...");

    DETRAY_VERBOSE_HOST(
        "Successfully built material maps for volume: " << this->name());

    // Construct the surfaces and give the volume to the next decorator
    return volume_decorator<detector_t>::build(det, ctx);
  }

 private:
  /// Check whether a surface with a given index @param sf_idx should receive
  /// material from this builder
  bool surface_has_map(const dindex sf_idx) const {
    return m_bin_data.contains(sf_idx);
  }

  /// Build the material grid for every surface that has material, add it to
  /// the detector @param det and set the global surface material link.
  ///
  /// @note The grids are built from the volume local masks in the builder,
  /// since the surfaces are not yet added to the detector
  void add_material_maps(detector_t& det) {
    DETRAY_VERBOSE_HOST("Build material maps for surfaces...");

    // The total number of surfaces that will be built by this builder
    const dindex n_surfaces{static_cast<dindex>(this->surfaces().size())};

    for (dindex sf_idx = 0u; sf_idx < n_surfaces; ++sf_idx) {
      if (!surface_has_map(sf_idx)) {
        continue;
      }

      auto& sf_desc = this->surfaces().at(sf_idx);

      DETRAY_DEBUG_HOST("-> surface #" << sf_idx << " sf_desc = " << sf_desc);

      // The axis spans of the material map (if empty, the extent of the
      // surface mask is used)
      darray<std::vector<scalar_type>, DIM> axis_spans{};
      if (auto axis_spans_itr = m_axis_spans.find(sf_idx);
          axis_spans_itr != m_axis_spans.end()) {
        axis_spans = axis_spans_itr->second;
      }

      // Construct the material map for a given surface shape and add it to
      // the detector (or find an identical map that is already there)
      auto [mat_id, mat_idx] =
          this->masks().template visit<detail::add_sf_material_map<material_t>>(
              sf_desc.mask(), m_factory, m_bin_data.at(sf_idx),
              m_n_bins.at(sf_idx), axis_spans, det._materials, m_deduplicate);

      // Make sure the material type was set correctly by the factory
      if (mat_id != sf_desc.material().id() || mat_idx == dindex_invalid) {
        std::stringstream err_msg;
        err_msg << "-> material id mismatch for surface " << sf_idx
                << ": expected " << sf_desc.material().id() << ", got "
                << mat_id;

        DETRAY_FATAL_HOST(err_msg.str());
        throw std::runtime_error(err_msg.str());
      }

      sf_desc.material().set_index(mat_idx);

      DETRAY_DEBUG_HOST("--> material link = " << sf_desc.material());
    }
  }

  /// Whether to share identical material maps between surfaces
  bool m_deduplicate{false};
  /// The surface this material map belongs to (index is volume local)
  std::map<dindex, std::vector<bin_data_type>> m_bin_data;
  /// Number of bins for the material grid axes
  std::map<dindex, darray<std::size_t, DIM>> m_n_bins{};
  /// The Axis spans for the material grid axes
  std::map<dindex, darray<std::vector<scalar_type>, DIM>> m_axis_spans{};
  /// Helper to generate empty grids
  mat_map_factory_t m_factory{};
};

namespace detail {

/// A functor to add a material map to a surface
///
/// If @c deduplicate is set, an identical material map that is already in the
/// material store is reused instead of adding a new one.
///
/// @returns the type id and index of the material map in the store
template <typename material_t>
struct add_sf_material_map {
  template <typename mask_coll_t, typename index_range_t,
            typename mat_factory_t, typename bin_data_t, std::size_t DIM,
            typename material_store_t, concepts::scalar scalar_t>
  DETRAY_HOST inline std::pair<typename material_t::id, dindex> operator()(
      [[maybe_unused]] const mask_coll_t& mask_coll,
      [[maybe_unused]] const index_range_t& index,
      [[maybe_unused]] const mat_factory_t& mat_factory,
      [[maybe_unused]] std::vector<bin_data_t>& bin_data,
      [[maybe_unused]] const darray<std::size_t, DIM>& n_bins,
      [[maybe_unused]] const darray<std::vector<scalar_t>, DIM>& axis_spans,
      [[maybe_unused]] material_store_t& mat_store,
      [[maybe_unused]] const bool deduplicate = false) const {
    using mask_t = typename mask_coll_t::value_type;

    // No material maps for line surfaces
    if constexpr (!concepts::line_object<mask_t> && mask_t::shape::dim == DIM) {
      // Map a grid onto the surface mask (the boundaries are taken from
      // the @c axis_spans variable, if it is not empty)
      mask_t sf_mask = {};
      if constexpr (concepts::interval<index_range_t>) {
        using index_t = typename index_range_t::index_type;

        // Find the true surface extent over all masks
        sf_mask = mask_coll.at(index.lower());

        if (index.size() > 1u) {
          const index_range_t other_masks{
              index.lower() + 1u, static_cast<index_t>(index.size() - 1u)};

          // Merge sub-masks
          for (const auto& sub_mask :
               detray::ranges::subrange(mask_coll, other_masks)) {
            sf_mask = sf_mask + sub_mask;
          }
        }
      } else {
        sf_mask = mask_coll.at(index);
      }

      auto mat_grid = mat_factory.new_grid(sf_mask, n_bins, {}, {}, axis_spans);

      // The detector only knows the non-owning grid types
      using non_owning_t = typename decltype(mat_grid)::template type<false>;

      // Not every mask shape might be used for material maps
      if constexpr (types::contains<material_t, non_owning_t>) {
        DETRAY_VERBOSE_HOST("Filling material grid...");

        // Add the material slabs to the grid
        for (const auto& bin : bin_data) {
          mat_grid.template populate<replace<>>(bin.local_bin_idx,
                                                bin.single_element);
        }

        constexpr auto gid{types::id<material_t, non_owning_t>};

        // Look for an identical material grid in the detector
        if (deduplicate) {
          const auto& grid_coll = mat_store.template get<gid>();
          for (dindex i = 0u; i < grid_coll.size(); ++i) {
            if (is_identical_grid(mat_grid, grid_coll[i])) {
              DETRAY_VERBOSE_HOST(
                  "Found identical material grid:" << gid << " at index " << i);
              return {gid, i};
            }
          }
        }

        // Add the material grid to the detector
        mat_store.template push_back<gid>(mat_grid);
        DETRAY_VERBOSE_HOST("Built material grid:" << gid << ":\n"
                                                   << mat_grid.axes());

        // Return the index of the new material map
        return {gid, static_cast<dindex>(mat_store.template size<gid>() - 1u)};
      } else {
        return {material_t::id::e_none, dindex_invalid};
      }
    } else {
      return {material_t::id::e_none, dindex_invalid};
    }
  }
};

}  // namespace detail

}  // namespace detray
