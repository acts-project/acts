// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/unmasked.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <memory>
#include <vector>

namespace detray {

/// @brief Factory class for material maps.
///
/// Uses a surface factory underneath the hood that handles the surface
/// construction. The material ID from the detector determines the type of
/// material map that will be produced. The factory is filled from [a vector] of
/// @c material_data .
///
/// @tparam detector_t type of detector that contains the material
/// @tparam index_t local bin index type
template <typename detector_t, typename index_t = dindex>
class material_map_factory final : public factory_decorator<detector_t> {
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;
  using index_type = index_t;

  using base_factory = factory_decorator<detector_t>;
  using placeholder_factory_t = surface_factory<detector_t, unmasked<>>;

 public:
  using scalar_type = dscalar<typename detector_t::algebra_type>;
  using data_type = material_data<scalar_type>;

  using base_factory::operator();

  /// Factory with surfaces potentially already filled or empty placeholder
  /// that will not be used.
  DETRAY_HOST
  explicit material_map_factory(
      std::unique_ptr<surface_factory_interface<detector_t>> sf_factory =
          std::make_unique<placeholder_factory_t>())
      : base_factory(std::move(sf_factory)) {}

  /// @returns the number of material instances that will be built by the
  /// factory
  DETRAY_HOST
  auto n_materials() const -> dindex {
    const dindex n_surfaces{static_cast<dindex>(m_links.size())};

    // Need exactly one material per surface
    assert(m_n_bins.size() == n_surfaces);
    assert(m_axis_spans.size() == n_surfaces);
    assert(m_materials.size() == n_surfaces);
    assert(m_thickness.size() == n_surfaces);

    return n_surfaces;
  }

  /// @returns the material links to the surfaces (counted for this volume)
  DETRAY_HOST
  auto links() const -> const
      std::map<std::size_t, std::pair<material_id, std::vector<index_type>>> & {
    return m_links;
  }

  /// @returns the raw materials that are currently in the factory
  DETRAY_HOST
  auto materials() const
      -> const std::map<std::size_t, std::vector<material<scalar_type>>> & {
    return m_materials;
  }

  /// @returns the material thickness currently held by the factory
  DETRAY_HOST
  auto thickness() const
      -> const std::map<std::size_t, std::vector<scalar_type>> & {
    return m_thickness;
  }

  /// Add all necessary components to the factory for a material map
  ///
  /// @param id the identifier for thr type of material map
  /// @param n_bins the number of bins for the material grid axes
  /// @param axis_spans the axis ranges for the material grid axes
  /// @param mat_data contains the material data
  /// @param indices local bin indices in the same order as the material data
  DETRAY_HOST
  void add_material(material_id id, data_type &&mat_data,
                    std::vector<std::size_t> &&n_bins,
                    std::vector<std::vector<scalar_type>> &&axis_spans,
                    std::vector<index_type> &&indices) {
    DETRAY_VERBOSE_HOST("Received material map data");

    auto [sf_index, mat, thickness] = std::move(mat_data).get_data();

    assert(!detail::is_invalid_value(sf_index));

    m_links.emplace(sf_index, std::make_pair(id, std::move(indices)));
    m_n_bins.emplace(sf_index, std::move(n_bins));
    m_axis_spans.emplace(sf_index, std::move(axis_spans));
    m_materials.emplace(sf_index, std::move(mat));
    m_thickness.emplace(sf_index, std::move(thickness));
  }

  /// Clear old data
  DETRAY_HOST
  auto clear() -> void override {
    m_links.clear();
    m_n_bins.clear();
    m_axis_spans.clear();
    m_materials.clear();
    m_thickness.clear();
  }

  /// @brief Add data for material maps to the containers of a volume builder.
  ///
  /// This assumes that the surfaces have already been added (e.g. by this
  /// factories underlying surface factory).
  ///
  /// @note This does not override the pure virtual function from the surface
  /// factory interface, but presents an overload for the case when material
  /// should be added.
  ///
  /// @param surfaces surface container of the volume builder that should get
  ///                 decorated with material maps.
  /// @param materials map the local bin indices and their content to a
  ///                  volume local surface index.
  template <typename bin_data_t, std::size_t N>
  DETRAY_HOST auto operator()(
      typename detector_t::surface_lookup_container &surfaces,
      std::map<dindex, std::vector<bin_data_t>> &material_map,
      std::map<dindex, darray<std::size_t, N>> &n_bins,
      std::map<dindex, darray<std::vector<scalar_type>, N>> &axis_spans) {
    DETRAY_VERBOSE_HOST("Add material maps...");

    using link_t = typename detector_t::surface_type::material_link;

    // Check data consistency
    const std::size_t n_materials{this->n_materials()};
    if (n_materials == 0) {
      return;
    }

    assert(surfaces.size() >= n_materials);

    DETRAY_VERBOSE_HOST("-> Have " << m_materials.size()
                                   << " configured material map(s) to assign");

    // Add the material only to those surfaces that the data links against
    for (auto &[i, materials] : m_materials) {
      const auto sf_idx{static_cast<dindex>(i)};

      DETRAY_DEBUG_HOST("--> #" << sf_idx << " sf = " << surfaces.at(sf_idx));

      // Copy the number of bins to the builder
      assert(m_n_bins.at(sf_idx).size() == N);
      n_bins.emplace(sf_idx, darray<std::size_t, N>{});
      std::ranges::copy_n(m_n_bins.at(sf_idx).begin(), N,
                          n_bins.at(sf_idx).begin());

      // Copy the axis spans to the builder (if present)
      axis_spans.emplace(sf_idx, darray<std::vector<scalar_type>, N>{});
      for (std::size_t in = 0; in < N; ++in) {
        if (m_axis_spans.at(sf_idx).size() > in) {
          axis_spans.at(sf_idx).at(in) = m_axis_spans.at(sf_idx).at(in);
        }
      }

      // Combine the material slab with its local bin index
      for (const auto [j, mat] : detray::views::enumerate(materials)) {
        scalar_type t{m_thickness.at(sf_idx)[j]};
        auto &loc_bin_idx{m_links.at(sf_idx).second.at(j)};

        material_slab<scalar_type> mat_slab{mat, t};
        bin_data_t data{loc_bin_idx, mat_slab};

        auto search = material_map.find(sf_idx);
        if (search == material_map.end()) {
          material_map.emplace(sf_idx, std::vector<bin_data_t>{data});
        } else {
          search->second.push_back(data);
        }
      }

      auto map_id{m_links.at(sf_idx).first};
      DETRAY_DEBUG_HOST("--> Set surface material id: " << map_id);
      surfaces.at(sf_idx).material() = link_t{map_id, dindex_invalid};
    }

    DETRAY_DEBUG_HOST("Finished");
  }

 private:
  /// Type and position(s) of the material in the detector material collection
  std::map<std::size_t, std::pair<material_id, std::vector<index_type>>>
      m_links{};
  /// Number of bins for the material grid axes
  std::map<std::size_t, std::vector<std::size_t>> m_n_bins{};
  /// Axis ranges for the material grid axes
  std::map<std::size_t, std::vector<std::vector<scalar_type>>> m_axis_spans{};
  /// Material thickness per surface
  std::map<std::size_t, std::vector<scalar_type>> m_thickness{};
  /// The pre-computed material to be wrapped in a slab per surface material
  /// map, sorted the same way as the bin index vector
  std::map<std::size_t, std::vector<material<scalar_type>>> m_materials{};
};

}  // namespace detray
