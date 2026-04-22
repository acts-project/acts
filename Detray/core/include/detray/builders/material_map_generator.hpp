// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/detail/surface_kernels.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/material/material.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/grid/axis.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace detray {

namespace detail {

/// Generate material along z bins for a cylinder material grid
template <concepts::scalar scalar_t>
inline std::vector<material_slab<scalar_t>> generate_cyl_mat(
    const dvector<dvector<scalar_t>> &bounds_vec, const std::size_t nbins,
    material<scalar_t> mat, const scalar_t t, const scalar_t scalor) {
  std::vector<material_slab<scalar_t>> ts;
  ts.reserve(nbins);

  // Find the surface bounds from the masks
  scalar_t z_lower{detail::invalid_value<scalar_t>()};
  scalar_t z_upper{-detail::invalid_value<scalar_t>()};
  for (const auto &bounds : bounds_vec) {
    z_lower = math::min(z_lower, bounds.at(cylinder2D::e_lower_z));
    z_upper = math::max(z_upper, bounds.at(cylinder2D::e_upper_z));
  }

  // Generate material steps (quadratic with distance to origin)
  const scalar_t length{math::fabs(z_upper - z_lower)};
  scalar_t z{z_lower};
  const scalar_t z_step{length / static_cast<scalar_t>(nbins - 1u)};
  for (std::size_t n = 0u; n < nbins; ++n) {
    ts.emplace_back(mat, static_cast<scalar_t>(scalor * z * z) + t);
    z += z_step;
  }

  return ts;
}

/// Generate material along r bins for a disc material grid
template <concepts::scalar scalar_t>
inline std::vector<material_slab<scalar_t>> generate_disc_mat(
    const dvector<dvector<scalar_t>> &bounds_vec, const std::size_t nbins,
    material<scalar_t> mat, const scalar_t t, const scalar_t scalor) {
  std::vector<material_slab<scalar_t>> ts;
  ts.reserve(nbins);

  // Find the surface bounds from the masks
  scalar_t r_inner{detail::invalid_value<scalar_t>()};
  scalar_t r_outer{0.f};
  for (const auto &bounds : bounds_vec) {
    r_inner = math::min(r_inner, bounds.at(ring2D::e_inner_r));
    r_outer = math::max(r_outer, bounds.at(ring2D::e_outer_r));
  }

  scalar_t r{r_inner};
  const scalar_t r_step{(r_outer - r_inner) /
                        static_cast<scalar_t>(nbins - 1u)};
  for (std::size_t n = 0u; n < nbins; ++n) {
    ts.emplace_back(mat, static_cast<scalar_t>(scalor * r) + t);
    r += r_step;
  }

  return ts;
}

}  // namespace detail

/// @brief Configuration for the material map generator
template <concepts::scalar scalar_t>
struct material_map_config {
  /// How to configure the generation for a specific type of material map
  struct map_config {
    /// Type id value of the material map in a given detector
    dindex map_id{dindex_invalid};
    /// Number of bins for material maps
    darray<std::size_t, 2> n_bins{20u, 20u};
    /// Along which of the two axes to scale the material thickness
    std::size_t axis_index{0u};
    /// Material to be filled into the maps
    material<scalar_t> mapped_material{silicon_tml<scalar_t>()};
    /// Minimal thickness of the material slabs in the material maps
    scalar_t thickness{0.15f * unit<scalar_t>::mm};
    /// Scale factor for the material thickness calculation
    scalar_t scalor{1.f};
    /// How to vary the material thickness along the bins
    std::function<std::vector<material_slab<scalar_t>>(
        const dvector<dvector<scalar_t>> &, const std::size_t,
        material<scalar_t>, const scalar_t, const scalar_t)>
        mat_generator{detray::detail::generate_cyl_mat<scalar_t>};
  };

  /// Configuration for different types of material maps:
  /// The key consists of the mask type id and the surface type id
  std::map<std::pair<unsigned int, surface_id>, map_config> m_map_configs;

  /// Setters
  /// @{
  template <typename mask_id_t>
  constexpr material_map_config &set_map_config(mask_id_t map_type,
                                                surface_id sf_type,
                                                const map_config &cfg) {
    const auto key =
        std::make_pair(static_cast<unsigned int>(map_type), sf_type);
    if (!m_map_configs.contains(key)) {
      m_map_configs.emplace(key, cfg);
    } else {
      m_map_configs.at(key) = cfg;
    }
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  template <typename mask_id_t>
  bool has_config(mask_id_t map_type, surface_id sf_type) const {
    const auto key =
        std::make_pair(static_cast<unsigned int>(map_type), sf_type);
    return m_map_configs.contains(key);
  }
  template <typename mask_id_t>
  const auto &get_map_config(mask_id_t map_type, surface_id sf_type) const {
    const auto key =
        std::make_pair(static_cast<unsigned int>(map_type), sf_type);
    return m_map_configs.at(key);
  }
  template <typename mask_id_t>
  auto &get_map_config(mask_id_t map_type, surface_id sf_type) {
    const auto key =
        std::make_pair(static_cast<unsigned int>(map_type), sf_type);
    return m_map_configs.at(key);
  }
  /// @}
};

/// @brief Surface factory decorator that adds material maps to surfaces
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t>
class material_map_generator final : public factory_decorator<detector_t> {
  using scalar_t = dscalar<typename detector_t::algebra_type>;

 public:
  /// Construct from configuration @param cfg
  DETRAY_HOST
  material_map_generator(
      std::unique_ptr<surface_factory_interface<detector_t>> factory,
      const material_map_config<scalar_t> cfg)
      : factory_decorator<detector_t>(std::move(factory)), m_cfg{cfg} {}

  /// Call the underlying surface factory and record the surface range that
  /// was produced
  ///
  /// @param volume the volume the portals need to be added to.
  /// @param surfaces the surface collection to wrap and to add the portals to
  /// @param transforms the transforms of the surfaces.
  /// @param masks the masks of the surfaces.
  /// @param ctx the geometry context (not needed for portals).
  DETRAY_HOST
  auto operator()(typename detector_t::volume_type &volume,
                  typename detector_t::surface_lookup_container &surfaces,
                  typename detector_t::transform_container &transforms,
                  typename detector_t::mask_container &masks,
                  typename detector_t::geometry_context ctx = {})
      -> dindex_range override {
    const dindex_range idx_range =
        (*this->get_factory())(volume, surfaces, transforms, masks, ctx);

    m_surface_range = idx_range;

    return idx_range;
  }

  /// Create material maps for all surfaces that the underlying surface
  /// factory builds.
  ///
  /// @param surfaces surface container of the volume builder that should get
  ///                 decorated with material maps.
  /// @param masks mask store of the volume builder that contains the shapes
  ///              and extents of the surfaces.
  /// @param material_map map the local bin indices and their content to a
  ///                  volume local surface index.
  /// @param n_bins the number of bins per axis per surface.
  template <typename bin_data_t, std::size_t N>
  DETRAY_HOST auto operator()(
      typename detector_t::surface_lookup_container &surfaces,
      const typename detector_t::mask_container &masks,
      std::map<dindex, std::vector<bin_data_t>> &material_map,
      std::map<dindex, darray<std::size_t, N>> &n_bins) {
    DETRAY_DEBUG_HOST("Generate material maps...");

    static_assert(N == 2u, "This generator only supports 2D material maps");

    using link_t = typename detector_t::surface_type::material_link;
    using mask_id = typename detector_t::masks::id;
    using material_id = typename detector_t::material::id;

    assert(surfaces.size() >= (m_surface_range[1] - m_surface_range[0]));

    // Add the material only to those surfaces that the data links against
    unsigned int i{0u};
    for (auto &sf : detray::ranges::subrange(surfaces, m_surface_range)) {
      const auto sf_idx{m_surface_range[0] + i};

      // Skip line surfaces, if any are defined
      if constexpr (types::contains<typename detector_t::material,
                                    material_rod<scalar_t>>) {
        const mask_id sf_mask_id = sf.mask().id();
        if (sf_mask_id == mask_id::e_straw_tube ||
            sf_mask_id == mask_id::e_drift_cell) {
          continue;
        }
      }

      // Get the correct configuration for this surface
      if (!m_cfg.has_config(sf.mask().id(), sf.id())) {
        continue;
      }
      const auto &map_cfg = m_cfg.get_map_config(sf.mask().id(), sf.id());

      // Copy the number of bins to the builder
      assert(map_cfg.n_bins.size() == N);
      n_bins.emplace(sf_idx, darray<std::size_t, N>{});
      std::ranges::copy_n(map_cfg.n_bins.begin(), N, n_bins.at(sf_idx).begin());

      // Scale material thickness either over e.g. r- or z-bins
      const std::size_t bins = map_cfg.n_bins.at(map_cfg.axis_index);

      using sf_kernels =
          detail::surface_kernels<typename detector_t::algebra_type>;

      const auto bounds_vec =
          masks.template visit<typename sf_kernels::get_mask_values>(sf.mask());

      // Generate the raw material
      const auto material{
          map_cfg.mat_generator(bounds_vec, bins, map_cfg.mapped_material,
                                map_cfg.thickness, map_cfg.scalor)};

      // Add the material slabs with their local bin indices to the
      // volume builder
      for (dindex bin0 = 0u; bin0 < map_cfg.n_bins.at(0); ++bin0) {
        for (dindex bin1 = 0u; bin1 < map_cfg.n_bins.at(1); ++bin1) {
          bin_data_t data{
              axis::multi_bin<N>{bin0, bin1},
              material.at((map_cfg.axis_index == 0u) ? bin0 : bin1)};

          auto search = material_map.find(sf_idx);
          if (search == material_map.end()) {
            material_map.emplace(sf_idx, std::vector<bin_data_t>{data});
          } else {
            search->second.push_back(data);
          }
        }
      }

      // Set the initial surface material link (will be updated when
      // added to the detector)
      sf.material() =
          link_t{static_cast<material_id>(map_cfg.map_id), dindex_invalid};
      i++;
    }
  }

 private:
  /// Material generator configuration
  material_map_config<scalar_t> m_cfg;
  /// Range of surface indices for which to generate material
  dindex_range m_surface_range{};
};

}  // namespace detray
