// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/concepts.hpp"
#include "detray/geometry/tracking_volume.hpp"

// System include(s)
#include <memory>
#include <string>
#include <string_view>

namespace detray {

template <concepts::detector detector_t>
class surface_factory_interface;

template <concepts::detector detector_t>
class volume_decorator;

/// Options for the volume builder classes
struct volume_builder_options {
  /// Globally deduplicate identical data that is added by a builder to the
  /// detector
  bool m_deduplicate{true};
  /// Globally sort data including the one that is added by a builder to the
  /// detector
  bool m_sort{true};

  /// Returns whether global deduplication is enabled.
  /// @returns The deduplication setting.
  /// @{
  constexpr bool deduplicate() const { return m_deduplicate; }
  /// Returns whether global sorting is enabled.
  /// @returns The sorting setting.
  constexpr bool sort() const { return m_sort; }
  /// @}

  /// Sets whether global deduplication is enabled.
  /// @param toggle New deduplication setting.
  /// @returns This options object.
  /// @{
  constexpr volume_builder_options &deduplicate(bool toggle) {
    m_deduplicate = toggle;
    return *this;
  }
  /// Sets whether global sorting is enabled.
  /// @param toggle New sorting setting.
  /// @returns This options object.
  constexpr volume_builder_options &sort(bool toggle) {
    m_sort = toggle;
    return *this;
  }
  /// @}
};

/// @brief Interface for volume builders (and volume builder decorators)
template <concepts::detector detector_t>
class volume_builder_interface {
  // Access protected methods
  friend class volume_decorator<detector_t>;

 public:
  /// Algebra type used by the detector.
  using algebra_type = typename detector_t::algebra_type;
  /// Scalar type used by the detector's algebra.
  using scalar_type = dscalar<algebra_type>;

  virtual ~volume_builder_interface() = default;

  /// @returns the global index for the volume
  /// @note the correct index is only available after calling @c init_vol
  DETRAY_HOST
  virtual auto vol_index() const -> dindex = 0;

  /// Toggles whether sensitive surfaces are added to the brute force method.
  /// @param toggle Whether to add sensitive surfaces to the brute force method.
  DETRAY_HOST
  virtual void has_accel(bool toggle) = 0;

  /// @returns Whether sensitive surfaces are added to the brute force method.
  DETRAY_HOST
  virtual bool has_accel() const = 0;

  /// Sets the name for the volume.
  /// @param volume_name Name to assign to the volume.
  DETRAY_HOST
  virtual void set_name(std::string volume_name) = 0;

  /// @returns the name of the volume
  DETRAY_HOST
  virtual std::string_view name() = 0;

  /// Returns const access to the volume being built.
  /// @returns The volume being built.
  DETRAY_HOST
  virtual auto operator()() const -> const
      typename detector_t::volume_type & = 0;
  /// Returns mutable access to the volume being built.
  /// @returns The volume being built.
  DETRAY_HOST
  virtual auto operator()() -> typename detector_t::volume_type & = 0;

  /// Adds a volume and all of its contents to a detector.
  /// @param det Detector receiving the volume.
  /// @param opt Options controlling deduplication and sorting.
  /// @param ctx Geometry context used while building.
  /// @returns Pointer to the built volume in the detector.
  DETRAY_HOST
  virtual auto build(detector_t &det, const volume_builder_options &opt,
                     typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * = 0;

  /// Adds a transform for the volume placement.
  /// @param trf Transform to copy into the placement.
  DETRAY_HOST
  virtual void add_volume_placement(
      const dtransform3D<algebra_type> &trf = {}) = 0;

  /// Adds a volume placement using a translation and identity rotation.
  /// @param t Translation of the volume placement.
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t) = 0;

  /// Adds a volume placement using a translation and its x- and z-axes.
  /// @param t Translation of the volume placement.
  /// @param x New x-axis of the placement.
  /// @param z New z-axis of the placement.
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t,
                                    const dvector3D<algebra_type> &x,
                                    const dvector3D<algebra_type> &z) = 0;

  /// Adds surfaces to the volume through a surface factory.
  /// @param sf_factory Factory that supplies surfaces for the volume.
  /// @param ctx Geometry context used while adding surfaces.
  DETRAY_HOST
  virtual void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) = 0;

 protected:
  /// Returns mutable access to the builder's surface lookup container.
  /// @returns The surface lookup container.
  /// @{
  virtual typename detector_t::surface_lookup_container &surfaces() = 0;
  /// Returns mutable access to the builder's transforms.
  /// @returns The transform container.
  virtual typename detector_t::transform_container &transforms() = 0;
  /// Returns mutable access to the builder's masks.
  /// @returns The mask container.
  virtual typename detector_t::mask_container &masks() = 0;
  /// @}
};

/// Decorator that forwards volume-building operations to an underlying builder.
///
/// Can be volume builders that introduce special sorting/memory layout, or
/// accelerator builders, like the grid builder.
template <concepts::detector detector_t>
class volume_decorator : public volume_builder_interface<detector_t> {
 public:
  /// Scalar type used by the detector's algebra.
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  /// Wraps a volume builder.
  /// @param vol_builder Builder to decorate; must not be null.
  DETRAY_HOST
  explicit volume_decorator(
      std::unique_ptr<volume_builder_interface<detector_t>> vol_builder)
      : m_builder(std::move(vol_builder)) {
    assert(m_builder != nullptr);
  }

  /// Returns the mutable volume being built.
  /// @returns The wrapped volume.
  DETRAY_HOST
  auto operator()() -> typename detector_t::volume_type & override {
    return m_builder->operator()();
  }

  /// Returns const access to the wrapped volume.
  /// @returns The wrapped volume.
  DETRAY_HOST
  auto operator()() const -> const typename detector_t::volume_type & override {
    return m_builder->operator()();
  }

  /// Returns the index assigned by the wrapped builder.
  /// @returns The global volume index.
  DETRAY_HOST
  auto vol_index() const -> dindex override { return m_builder->vol_index(); }

  /// Forwards the brute-force acceleration setting to the wrapped builder.
  /// @param toggle Whether sensitive surfaces are added to brute force lookup.
  DETRAY_HOST
  void has_accel(bool toggle) override { m_builder->has_accel(toggle); };

  /// Returns whether sensitive surfaces are added to brute force lookup.
  /// @returns The brute-force lookup setting.
  DETRAY_HOST
  bool has_accel() const override { return m_builder->has_accel(); }

  /// Sets the wrapped volume's name.
  /// @param volume_name Name to assign.
  DETRAY_HOST
  void set_name(std::string volume_name) override {
    return m_builder->set_name(volume_name);
  }

  /// Returns the wrapped volume's name.
  /// @returns The volume name.
  DETRAY_HOST std::string_view name() override { return m_builder->name(); }

  /// Builds the wrapped volume in the detector.
  /// @param det Detector receiving the volume.
  /// @param opt Options controlling deduplication and sorting.
  /// @returns Pointer to the built volume.
  DETRAY_HOST
  auto build(detector_t &det, const volume_builder_options &opt,
             typename detector_t::geometry_context /*ctx*/ = {}) ->
      typename detector_t::volume_type * override {
    return m_builder->build(det, opt);
  }

  /// Adds a placement transform to the wrapped builder.
  /// @param trf Transform to copy into the placement.
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::transform3_type &trf = {}) override {
    return m_builder->add_volume_placement(trf);
  }

  /// Adds a placement using a translation and identity rotation.
  /// @param t Translation of the placement.
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type &t) override {
    return m_builder->add_volume_placement(t);
  }

  /// Adds a placement using a translation and its x- and z-axes.
  /// @param t Translation of the placement.
  /// @param x New x-axis of the placement.
  /// @param z New z-axis of the placement.
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type &t,
      const typename detector_t::vector3_type &x,
      const typename detector_t::vector3_type &z) override {
    return m_builder->add_volume_placement(t, x, z);
  }

  /// Adds surfaces through the wrapped builder's factory interface.
  /// @param sf_factory Factory that supplies surfaces.
  /// @param ctx Geometry context forwarded to the builder.
  DETRAY_HOST
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) override {
    return m_builder->add_surfaces(std::move(sf_factory), ctx);
  }

 protected:
  /// Access to the underlying builder and its data.
  /// @{
  /// Returns a pointer to the wrapped builder.
  /// @returns Pointer to the wrapped builder.
  volume_builder_interface<detector_t> *get_builder() {
    return m_builder.get();
  }
  /// Returns mutable access to the wrapped builder's surface lookup container.
  /// @returns The surface lookup container.
  typename detector_t::surface_lookup_container &surfaces() override {
    return m_builder->surfaces();
  }
  /// Returns mutable access to the wrapped builder's transforms.
  /// @returns The transform container.
  typename detector_t::transform_container &transforms() override {
    return m_builder->transforms();
  }
  /// Returns mutable access to the wrapped builder's masks.
  /// @returns The mask container.
  typename detector_t::mask_container &masks() override {
    return m_builder->masks();
  }
  /// @}

 private:
  std::unique_ptr<volume_builder_interface<detector_t>> m_builder;
};

}  // namespace detray
