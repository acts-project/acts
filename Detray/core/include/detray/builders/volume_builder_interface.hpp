// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/geometry/tracking_volume.hpp"

// System include(s)
#include <memory>
#include <string>
#include <string_view>

namespace detray {

template <typename detector_t>
class surface_factory_interface;

template <typename detector_t>
class volume_decorator;

/// @brief Interface for volume builders (and volume builder decorators)
template <typename detector_t>
class volume_builder_interface {
  // Access protected methods
  friend class volume_decorator<detector_t>;

 public:
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;

  virtual ~volume_builder_interface() = default;

  /// @returns the global index for the volume
  /// @note the correct index is only available after calling @c init_vol
  DETRAY_HOST
  virtual auto vol_index() const -> dindex = 0;

  /// Toggles whether sensitive surfaces are added to the brute force method
  DETRAY_HOST
  virtual void has_accel(bool toggle) = 0;

  /// @returns whether sensitive surfaces are added to the brute force method
  DETRAY_HOST
  virtual bool has_accel() const = 0;

  /// Sets the name @param volume_name for the volume
  DETRAY_HOST
  virtual void set_name(std::string volume_name) = 0;

  /// @returns the name of the volume
  DETRAY_HOST
  virtual std::string_view name() = 0;

  /// @returns reading access to the volume
  DETRAY_HOST
  virtual auto operator()() const -> const
      typename detector_t::volume_type & = 0;
  DETRAY_HOST
  virtual auto operator()() -> typename detector_t::volume_type & = 0;

  /// @brief Adds a volume and all of its contents to a detector
  DETRAY_HOST
  virtual auto build(detector_t &det,
                     typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * = 0;

  /// @brief Add the transform for the volume placement - copy
  DETRAY_HOST
  virtual void add_volume_placement(
      const dtransform3D<algebra_type> &trf = {}) = 0;

  /// @brief Add the transform for the volume placement from the translation
  /// @param t. The rotation will be the identity matrix.
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t) = 0;

  /// @brief Add the transform for the volume placement from the translation
  /// @param t , the new x- and z-axis (@param x, @param z).
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t,
                                    const dvector3D<algebra_type> &x,
                                    const dvector3D<algebra_type> &z) = 0;

  /// @brief Add surfaces to the volume
  /// @returns the index range of the sensitives in the temporary surface
  /// container used by the factory ( gets final update in @c build() )
  DETRAY_HOST
  virtual void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) = 0;

 protected:
  /// Access to builder data
  /// @{
  virtual typename detector_t::surface_lookup_container &surfaces() = 0;
  virtual typename detector_t::transform_container &transforms() = 0;
  virtual typename detector_t::mask_container &masks() = 0;
  /// @}
};

/// @brief Decorator for the volume builder.
///
/// Can be volume builders that introduce special sorting/memory layout, or
/// accelerator builders, like the grid builder.
template <typename detector_t>
class volume_decorator : public volume_builder_interface<detector_t> {
 public:
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  DETRAY_HOST
  explicit volume_decorator(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : m_builder(std::move(vol_builder)) {
    assert(m_builder != nullptr);
  }

  DETRAY_HOST
  auto operator()() -> typename detector_t::volume_type & override {
    return m_builder->operator()();
  }

  DETRAY_HOST
  auto operator()() const -> const typename detector_t::volume_type & override {
    return m_builder->operator()();
  }

  DETRAY_HOST
  auto vol_index() const -> dindex override { return m_builder->vol_index(); }

  DETRAY_HOST
  void has_accel(bool toggle) override { m_builder->has_accel(toggle); };

  DETRAY_HOST
  bool has_accel() const override { return m_builder->has_accel(); }

  DETRAY_HOST
  void set_name(std::string volume_name) override {
    return m_builder->set_name(volume_name);
  }

  DETRAY_HOST std::string_view name() override { return m_builder->name(); }

  DETRAY_HOST
  auto build(detector_t &det,
             typename detector_t::geometry_context /*ctx*/ = {}) ->
      typename detector_t::volume_type * override {
    return m_builder->build(det);
  }

  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::transform3_type &trf = {}) override {
    return m_builder->add_volume_placement(trf);
  }

  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type &t) override {
    return m_builder->add_volume_placement(t);
  }

  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type &t,
      const typename detector_t::vector3_type &x,
      const typename detector_t::vector3_type &z) override {
    return m_builder->add_volume_placement(t, x, z);
  }

  DETRAY_HOST
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) override {
    return m_builder->add_surfaces(std::move(sf_factory), ctx);
  }
  /// @}

 protected:
  /// Access to underlying builder
  /// @{
  volume_builder_interface<detector_t> *get_builder() {
    return m_builder.get();
  }
  typename detector_t::surface_lookup_container &surfaces() override {
    return m_builder->surfaces();
  }
  typename detector_t::transform_container &transforms() override {
    return m_builder->transforms();
  }
  typename detector_t::mask_container &masks() override {
    return m_builder->masks();
  }
  /// @}

 private:
  std::unique_ptr<volume_builder_interface<detector_t>> m_builder;
};

}  // namespace detray
