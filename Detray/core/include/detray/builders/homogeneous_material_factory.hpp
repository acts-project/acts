// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/surface_factory.hpp"
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/unmasked.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace detray {

/// @brief Bind components for material together.
template <concepts::scalar scalar_t>
class material_data {
 public:
  /// Construct empty data for a given surface
  ///
  /// @param sf_idx the index of the surface this material belongs to, needs
  ///               to be passed only if a special oredering must be observed
  DETRAY_HOST
  explicit constexpr material_data(
      const std::size_t sf_idx = detail::invalid_value<std::size_t>())
      : m_sf_index{sf_idx} {}

  /// Construct from a predefined material
  ///
  /// @param mat predefined material, see
  ///            'detray/material/predefined_materials.hpp'
  /// @param thickness of the material slab/rod
  /// @param sf_idx the index of the surface this material belongs to, needs
  ///               to be passed only if a special oredering must be observed
  DETRAY_HOST
  constexpr material_data(
      const scalar_t thickness, const material<scalar_t> &mat,
      const std::size_t sf_idx = detail::invalid_value<std::size_t>())
      : m_sf_index{sf_idx}, m_mat{mat}, m_thickness{thickness} {}

  /// Construct from all parameters:
  ///
  /// @param thickness of the material slab/rod
  /// @param material_paramters0 x0 is the radiation length
  /// @param material_paramters1 l0 is the nuclear interaction length
  /// @param material_paramters2 ar is the relative atomic mass
  /// @param material_paramters3 z is the nuclear charge number
  /// @param material_paramters4 molarRho is the molar density
  /// @param state of the material (liquid, solid etc.)
  /// @param sf_idx the index of the surface this material belongs to, needs
  ///               to be passed only if a special oredering must be observed
  DETRAY_HOST
  constexpr material_data(
      const scalar_t thickness, const std::vector<scalar_t> &material_paramters,
      const material_state state = material_state::e_solid,
      const std::size_t sf_idx = detail::invalid_value<std::size_t>())
      : m_sf_index{sf_idx},
        m_mat{material<scalar_t>{
            material_paramters.at(0), material_paramters.at(1),
            material_paramters.at(2), material_paramters.at(3),
            material_paramters.at(4), state}},
        m_thickness{thickness} {}

  /// @returns tuple based access to the contained material data.
  DETRAY_HOST
  constexpr auto get_data()
      -> std::tuple<std::size_t &, std::vector<material<scalar_t>> &,
                    std::vector<scalar_t> &> {
    return std::tie(m_sf_index, m_mat, m_thickness);
  }

  /// Append new material
  ///
  /// @param thickness of the material slab/rod
  /// @param mat predefined material, see
  ///            'detray/material/predefined_materials.hpp'
  DETRAY_HOST
  void append(const scalar_t thickness, const material<scalar_t> &mat) {
    m_mat.push_back(mat);
    m_thickness.push_back(thickness);
  }

  /// Append new material from @param other @c material_data - move
  DETRAY_HOST
  void append(material_data &&other) {
    m_mat.reserve(m_mat.size() + other.m_mat.size());
    std::ranges::move(other.m_mat, std::back_inserter(m_mat));

    m_thickness.reserve(m_thickness.size() + other.m_thickness.size());
    std::ranges::move(other.m_thickness, std::back_inserter(m_thickness));
  }

  /// Append new material
  ///
  /// @param thickness of the material slab/rod
  /// @param material_paramters0 x0 is the radiation length
  /// @param material_paramters1 l0 is the nuclear interaction length
  /// @param material_paramters2 ar is the relative atomic mass
  /// @param material_paramters3 z is the nuclear charge number
  /// @param material_paramters4 molarRho is the molar density
  /// @param state of the material (liquid, solid etc.)
  DETRAY_HOST
  void append(const scalar_t thickness,
              const std::vector<scalar_t> &material_paramters,
              const material_state state = material_state::e_solid) {
    m_mat.push_back(
        material<scalar_t>{material_paramters.at(0), material_paramters.at(1),
                           material_paramters.at(2), material_paramters.at(3),
                           material_paramters.at(4), state});
    m_thickness.push_back(thickness);
  }

  /// Output stream operator for material_data
  DETRAY_HOST friend std::ostream &operator<<(std::ostream &os,
                                              const material_data &mat_data) {
    os << "material_data{\nsf_index: ";
    if (detail::is_invalid_value(mat_data.m_sf_index)) {
      os << "invalid";
    } else {
      os << mat_data.m_sf_index;
    }
    os << ", materials: [\n";
    for (std::size_t i = 0; i < mat_data.m_mat.size(); ++i) {
      if (i > 0) {
        os << ",\n";
      }
      os << mat_data.m_mat.at(i);
    }
    os << "],\nthickness: [\n";
    for (std::size_t i = 0; i < mat_data.m_thickness.size(); ++i) {
      if (i > 0) {
        os << ", ";
      }
      os << mat_data.m_thickness.at(i);
    }
    os << "]}";
    return os;
  }

 private:
  /// Volume local index of the surface this material belongs to
  std::size_t m_sf_index{detail::invalid_value<std::size_t>()};
  /// The material parametrization
  std::vector<material<scalar_t>> m_mat{};
  /// Thickness/radius of the material slab/rod
  std::vector<scalar_t> m_thickness{};
};

/// @brief Factory class for homogeneous material.
///
/// Uses a surface factory underneath the hood that handles the surface
/// construction. The material ID from the detector determines the type of
/// material that will be produced. The factory is filled from [a vector] of
/// @c material_data .
///
/// @tparam detector_t type of detector that contains the material
template <typename detector_t>
class homogeneous_material_factory final
    : public factory_decorator<detector_t> {
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  using base_factory = factory_decorator<detector_t>;
  using placeholder_factory_t = surface_factory<detector_t, unmasked<>>;

 public:
  using scalar_type = dscalar<typename detector_t::algebra_type>;

  using base_factory::operator();

  /// Factory with surfaces potentially already filled or empty placeholder
  /// that will not be used.
  DETRAY_HOST
  explicit homogeneous_material_factory(
      std::unique_ptr<surface_factory_interface<detector_t>> sf_factory =
          std::make_unique<placeholder_factory_t>())
      : base_factory(std::move(sf_factory)) {
    static_assert(concepts::has_material_slabs<detector_t> ||
                      concepts::has_material_rods<detector_t>,
                  "No homogeneous surface material in detector type");
  }

  /// @returns the number of material instances that will be built by the
  /// factory
  DETRAY_HOST
  auto n_materials() const -> dindex {
    const dindex n_surfaces{static_cast<dindex>(m_links.size())};

    // Need exactly one material per surface
    assert(m_indices.empty() || (m_indices.size() == n_surfaces));
    assert(m_links.size() == n_surfaces);
    assert(m_materials.size() == n_surfaces);
    assert(m_thickness.size() == n_surfaces);

    return n_surfaces;
  }

  /// @returns the material links to the surfaces (counted for this volume)
  DETRAY_HOST
  auto links() const -> const std::vector<material_id> & { return m_links; }

  /// @returns the raw materials that are currently in the factory
  DETRAY_HOST
  auto materials() const -> const std::vector<material<scalar_type>> & {
    return m_materials;
  }

  /// @returns the material thickness currently held by the factory
  DETRAY_HOST
  auto thickness() const -> const std::vector<scalar_type> & {
    return m_thickness;
  }

  /// Add all necessary components to the factory for a single material slab
  /// or rod (determined by the @param id)
  DETRAY_HOST
  void add_material(material_id id, material_data<scalar_type> &&mat_data,
                    std::size_t index = detail::invalid_value<std::size_t>()) {
    DETRAY_DEBUG_HOST("homogeneous_material_factory::add_material single");

    auto [sf_index, mat, thickness] = std::move(mat_data).get_data();

    DETRAY_DEBUG_HOST("- sf_index=" << sf_index);
    DETRAY_DEBUG_HOST("- mat=" << mat.at(0));
    DETRAY_DEBUG_HOST("- thickness=" << thickness.at(0));

    // Only one homogeneous material slab/rod per surface
    assert(mat.size() == 1u);
    assert(thickness.size() == 1u);

    m_links.push_back(std::make_pair(id, static_cast<dindex>(index)));
    m_indices.push_back(sf_index);
    m_materials.push_back(mat.at(0));
    m_thickness.push_back(thickness.at(0));
  }

  /// Add all necessary components to the factory for multiple material slabs
  /// or rods (determined by the @param id)
  DETRAY_HOST
  void add_material(material_id id,
                    std::vector<material_data<scalar_type>> &&mat_data_vec) {
    DETRAY_DEBUG_HOST("homogeneous_material_factory::add_material multiple");
    // Read the material containers
    m_links.reserve(m_links.size() + mat_data_vec.size());
    m_materials.reserve(m_materials.size() + mat_data_vec.size());
    m_thickness.reserve(m_thickness.size() + mat_data_vec.size());

    // Add the material components
    for (auto &mat_data : mat_data_vec) {
      this->add_material(id, std::move(mat_data));
    }
  }

  /// Clear old data
  DETRAY_HOST
  auto clear() -> void override {
    m_links.clear();
    m_materials.clear();
    m_thickness.clear();
  }

  /// @brief Add material to the containers of a volume builder.
  ///
  /// This assumes that the surfaces have already been added (e.g. by this
  /// factories underlying surface factory). The material from this factory is
  /// added to the corresponding number of surfaces at the end of the calling
  /// volume builder.
  ///
  /// @note This does not override the pure virtual function from the surface
  /// factory interface, but presents an overload for the case when material
  /// should be added.
  ///
  /// @param surfaces surface container of the volume builder that should get
  ///                 decorated with material.
  /// @param material material store of the volume builder that the new
  ///                 materials get added to.
  DETRAY_HOST
  auto operator()(typename detector_t::surface_lookup_container &surfaces,
                  typename detector_t::material_container &materials) {
    DETRAY_VERBOSE_HOST("Add homogeneous material...");

    using link_t = typename detector_t::surface_type::material_link;
    DETRAY_DEBUG_HOST("-> Accumulated " << m_materials.size() << " materials");

    if (m_materials.empty()) {
      return;
    }

    DETRAY_DEBUG_HOST("-> n_materials=" << n_materials());
    DETRAY_DEBUG_HOST("-> links:");
    for ([[maybe_unused]] const auto &[i, link] :
         detray::views::enumerate(m_links)) {
      DETRAY_DEBUG_HOST("-> #" << i << ": " << link.first << " "
                               << link.second);
    }

    // Check that the surfaces were set up correctly
    const std::size_t n_materials{this->n_materials()};
    assert(surfaces.size() >= n_materials);

    // If no concrete surface ordering was passed, use index sequence
    // and add the materials to the trailing elements in the surfaces cont.
    if (m_indices.empty() ||
        std::ranges::find(m_indices, detail::invalid_value<std::size_t>()) !=
            m_indices.end()) {
      DETRAY_DEBUG_HOST("Indices were empty or invalid, recalculating");
      m_indices.resize(n_materials);
      std::iota(std::begin(m_indices), std::end(m_indices),
                surfaces.size() - n_materials);
    }
    DETRAY_DEBUG_HOST("-> Indices: " << DETRAY_LOG_VECTOR(m_indices));

    // Correctly index the data in this factory
    std::size_t sf_offset{*std::ranges::min_element(m_indices)};

    DETRAY_DEBUG_HOST("-> sf_offset=" << sf_offset);

    // Add the material to the surfaces that the data links against
    for (auto [i, sf] : detray::views::pick(surfaces, m_indices)) {
      if (sf.material().id() < detector_t::material::id::e_none) {
        DETRAY_WARN_HOST("Surface descriptor already has a material link:"
                         << sf.material().id()
                         << " skipping configured homogeneous assignment");
        continue;
      }

      std::size_t sf_idx{i - sf_offset};
      const material<scalar_type> &mat = m_materials.at(sf_idx);
      scalar_type t = m_thickness.at(sf_idx);

      DETRAY_DEBUG_HOST("-> Adding: sf_idx=" << sf_idx << ", i=" << i
                                             << ", sf=" << sf);
      DETRAY_DEBUG_HOST("           mat=" << mat << " thickness=" << t);

      dindex mat_idx{0u};

      if constexpr (concepts::has_material_slabs<detector_t>) {
        if (m_links.at(sf_idx).first == material_id::e_material_slab) {
          auto &mat_coll =
              materials.template get<material_id::e_material_slab>();

          material_slab<scalar_type> mat_slab{mat, t};
          mat_idx = this->insert_in_container(mat_coll, mat_slab,
                                              m_links.at(sf_idx).second);
        }
      }
      if constexpr (concepts::has_material_rods<detector_t>) {
        if (m_links.at(sf_idx).first == material_id::e_material_rod) {
          auto &mat_coll =
              materials.template get<material_id::e_material_rod>();

          material_rod<scalar_type> mat_rod{mat, t};
          mat_idx = this->insert_in_container(mat_coll, mat_rod,
                                              m_links.at(sf_idx).second);
        }
      }
      DETRAY_DEBUG_HOST("-> After insert: mat_idx=" << mat_idx);
      link_t new_link{m_links.at(sf_idx).first, mat_idx};

      auto &sf_desc = surfaces.at(static_cast<dindex>(i));
      DETRAY_DEBUG_HOST("-> Existing link: " << sf_desc.material());
      DETRAY_DEBUG_HOST("-> Setting new link: " << new_link);

      // Set the initial surface material link (will be updated when
      // added to the detector)
      sf_desc.material() = new_link;

      DETRAY_DEBUG_HOST("-> Added material to surface: " << sf);
    }
  }

 private:
  /// Material links of surfaces
  std::vector<std::pair<material_id, dindex>> m_links{};
  /// Position of the material in the detector material collection
  std::vector<std::size_t> m_indices{};
  /// Material thickness
  std::vector<scalar_type> m_thickness{};
  /// The pre-computed material to be wrapped in a slab or rod
  std::vector<material<scalar_type>> m_materials{};
};

}  // namespace detray
