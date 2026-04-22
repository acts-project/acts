// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: Remove this when gcc fixes their false positives.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wnull-dereference"
#endif

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace detray {

/// @brief Bind transform and mask data together for surface building.
///
/// Surface data is kept in separate detector containers and linked by indices,
/// so during geometry building, great care has to be taken to make sure that
/// all components of a surface get sorted and linked into the containers
/// together and associated with the correct surface.
template <typename detector_t>
class surface_data {
 public:
  using boundary_coll_type = std::vector<typename detector_t::scalar_type>;
  using navigation_link = typename detector_t::surface_type::navigation_link;

  /// Parametrized constructor - surface with only one mask
  ///
  /// @param type the surface type (portal|sensitive|passive)
  /// @param trf the surface placement transformation
  /// @param volume_link the mask volume link (used for navigation)
  /// @param mask_boundaries define the extent of the surface
  /// @param idx the index of the surface in the global detector lookup, needs
  ///            to be passed only if a special ordering should be observed
  DETRAY_HOST
  surface_data(
      const surface_id type, const typename detector_t::transform3_type &trf,
      navigation_link volume_link, const boundary_coll_type &mask_boundaries,
      const dindex idx = dindex_invalid,
      const std::uint64_t source = detail::invalid_value<std::uint64_t>())
      : m_type{type},
        m_volume_link{std::vector{volume_link}},
        m_index{idx},
        m_source{source},
        m_boundaries{std::vector{mask_boundaries}},
        m_transform{trf} {}

  /// Parametrized constructor
  ///
  /// @param type the surface type (portal|sensitive|passive)
  /// @param trf the surface placement transformation
  /// @param volume_link the volume link of each mask (used for navigation)
  /// @param mask_boundaries define the extent of each mask of the surface
  /// @param idx the index of the surface in the global detector lookup, needs
  ///            to be passed only if a special ordering should be observed
  DETRAY_HOST
  surface_data(
      const surface_id type, const typename detector_t::transform3_type &trf,
      std::vector<navigation_link> volume_link,
      const std::vector<boundary_coll_type> &mask_boundaries,
      const dindex idx = dindex_invalid,
      const std::uint64_t source = detail::invalid_value<std::uint64_t>())
      : m_type{type},
        m_volume_link{volume_link},
        m_index{idx},
        m_source{source},
        m_boundaries{mask_boundaries},
        m_transform{trf} {}

  /// Access the contained data through structured binding
  DETRAY_HOST
  auto get_data()
      -> std::tuple<surface_id &, std::vector<navigation_link> &, dindex &,
                    std::uint64_t &, std::vector<boundary_coll_type> &,
                    typename detector_t::transform3_type &> {
    return std::tie(m_type, m_volume_link, m_index, m_source, m_boundaries,
                    m_transform);
  }

 private:
  /// Surface type
  surface_id m_type;
  /// The index of the volume that this surface links to (per mask)
  std::vector<navigation_link> m_volume_link;
  /// The position of the surface in the detector containers, used to match
  /// the surface to e.g. its material
  dindex m_index;
  /// Source link (ACTS geoID)
  std::uint64_t m_source;
  /// Vector of mask boundary value collections
  std::vector<boundary_coll_type> m_boundaries;
  /// The surface placement
  typename detector_t::transform3_type m_transform;
};

/// @brief How to generate surfaces with their corresponding masks and
/// transforms.
///
/// Can be hard coded surface generation or json reader.
template <typename detector_t>
class surface_factory_interface {
 public:
  using navigation_link = typename detector_t::surface_type::navigation_link;

  surface_factory_interface() = default;
  virtual ~surface_factory_interface() = default;

  /// @returns the number of surfaces the factory will produce
  DETRAY_HOST
  virtual dindex size() const = 0;

  /// Add data to the factory
  /// @{
  DETRAY_HOST
  virtual void push_back(surface_data<detector_t> &&) = 0;

  DETRAY_HOST
  virtual void push_back(std::vector<surface_data<detector_t>> &&) = 0;
  /// @}

  /// Clear all data in the factory
  DETRAY_HOST
  virtual void clear() = 0;

  /// Construct detector components from the data in the factory and add them
  /// the containers of a volume builder.
  DETRAY_HOST
  virtual auto operator()(
      typename detector_t::volume_type &volume,
      typename detector_t::surface_lookup_container &surfaces,
      typename detector_t::transform_container &transforms,
      typename detector_t::mask_container &masks,
      typename detector_t::geometry_context ctx = {}) -> dindex_range = 0;

 protected:
  /// Insert a value in a container at a specific index
  ///
  /// @param cont the container to be filled
  /// @param value the value to be put into the container
  /// @param idx the position where to insert the value
  /// @param args optional additional parameters for the container access
  ///
  /// @returns the position where the value has been copied.
  template <typename container_t, typename... Args>
  DETRAY_HOST dindex insert_in_container(
      container_t &cont, const typename container_t::value_type value,
      const dindex idx, Args &&...args) const {
    // If no valid position is given, perform push back
    if (detail::is_invalid_value(idx)) {
      cont.push_back(value, std::forward<Args>(args)...);

      return static_cast<dindex>(cont.size() - 1u);
    } else {
      // Make sure the container size encompasses the new value
      if (cont.size() < idx + 1) {
        cont.resize(idx + 1, std::forward<Args>(args)...);
      }

      cont.at(idx, std::forward<Args>(args)...) = value;

      return idx;
    }
  }
};

/// @brief Decorator for the surface factories.
///
/// Delegates all mandatory method calls to the underlying surface factory
template <typename detector_t>
class factory_decorator : public surface_factory_interface<detector_t> {
 public:
  DETRAY_HOST
  explicit factory_decorator(
      std::unique_ptr<surface_factory_interface<detector_t>> factory)
      : m_factory(std::move(factory)) {
    if (m_factory == nullptr || m_factory.get() == nullptr) {
      throw std::runtime_error(
          "Surface factory decorator constructed with invalid base "
          "factory");
    }
  }

  /// @returns access to the underlying factory - const
  DETRAY_HOST
  const surface_factory_interface<detector_t> *get_factory() const {
    return m_factory.get();
  }

  /// Overwrite interface functions using callbacks to the base factory
  /// @{
  DETRAY_HOST
  dindex size() const override { return m_factory->size(); }

  DETRAY_HOST
  void push_back(surface_data<detector_t> &&data) override {
    m_factory->push_back(std::move(data));
  }

  DETRAY_HOST
  void push_back(std::vector<surface_data<detector_t>> &&data) override {
    m_factory->push_back(std::move(data));
  }

  DETRAY_HOST
  void clear() override { m_factory->clear(); }

  DETRAY_HOST
  auto operator()(typename detector_t::volume_type &volume,
                  typename detector_t::surface_lookup_container &surfaces,
                  typename detector_t::transform_container &transforms,
                  typename detector_t::mask_container &masks,
                  typename detector_t::geometry_context ctx = {})
      -> dindex_range override {
    return (*m_factory)(volume, surfaces, transforms, masks, ctx);
  }
  /// @}

 protected:
  /// @returns access to the underlying factory - non-const
  DETRAY_HOST
  surface_factory_interface<detector_t> *get_factory() {
    return m_factory.get();
  }

 private:
  /// Wrapped surface factory that new functionalit is added to by decoration
  std::unique_ptr<surface_factory_interface<detector_t>> m_factory;
};

}  // namespace detray
