// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/builders/deduplication_cache.hpp"
#include "detray/core/concepts.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace detray {

template <concepts::detector detector_t>
class surface_factory_interface;

template <concepts::detector detector_t>
class volume_decorator;

/// Options for the volume builder classes
struct volume_builder_options {
  /// Globally deduplicate identical data that is added by a builder to the
  /// detector (transforms, masks and material)
  /// @note Surfaces that share a transform cannot be aligned independently.
  /// The builder warns if this happens for sensitive surfaces.
  bool m_deduplicate{true};
  /// Globally sort data including the one that is added by a builder to the
  /// detector
  bool m_sort{true};

  /// Getters
  /// @returns The deduplication setting.
  /// @{
  constexpr bool deduplicate() const { return m_deduplicate; }
  /// Returns whether global sorting is enabled.
  /// @returns The sorting setting.
  constexpr bool sort() const { return m_sort; }
  /// @}

  /// Setters
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

  /// Print the volume builder configuration
  friend std::ostream &operator<<(std::ostream &out,
                                  const volume_builder_options &opt) {
    out << std::boolalpha;
    out << "\nBuild Options\n"
        << "----------------------------\n"
        << "  Deduplicate           : " << opt.deduplicate() << "\n"
        << "  Sort                  : " << opt.sort() << "\n";
    out << std::noboolalpha;

    return out;
  }
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

  /// Toggles whether sensitive surfaces are added to the brute force method
  /// @param toggle Whether to add sensitive surfaces to the brute force method.
  DETRAY_HOST
  virtual void has_accel(bool toggle) = 0;

  /// @returns whether sensitive surfaces are added to the brute force method
  DETRAY_HOST
  virtual bool has_accel() const = 0;

  /// Sets the name @p volume_name for the volume
  /// @param volume_name Name to assign to the volume.
  DETRAY_HOST
  virtual void set_name(std::string volume_name) = 0;

  /// @returns the name of the volume
  DETRAY_HOST
  virtual std::string_view name() = 0;

  /// Returns const access to the volume being built.
  /// @returns reading access to the volume
  DETRAY_HOST
  virtual auto operator()() const -> const
      typename detector_t::volume_type & = 0;
  /// Returns mutable access to the volume being built.
  /// @returns The volume being built.
  DETRAY_HOST
  virtual auto operator()() -> typename detector_t::volume_type & = 0;

  /// @brief Adds a volume and all of its contents to a detector
  /// @param det Detector receiving the volume.
  /// @param opt Options controlling deduplication and sorting.
  /// @param ctx Geometry context used while building.
  /// @returns Pointer to the built volume in the detector.
  DETRAY_HOST
  virtual auto build(detector_t &det, const volume_builder_options &opt,
                     typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * = 0;

  /// @brief Add the transform for the volume placement - copy
  /// @param trf Transform to copy into the placement.
  DETRAY_HOST
  virtual void add_volume_placement(
      const dtransform3D<algebra_type> &trf = {}) = 0;

  /// @brief Add the transform for the volume placement from the translation
  /// @p t. The rotation will be the identity matrix.
  /// @param t Translation of the volume placement.
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t) = 0;

  /// @brief Add the transform for the volume placement from the translation
  /// @p t , the new x- and z-axis (@p x, @p z).
  /// @param t Translation of the volume placement.
  /// @param x New x-axis of the placement.
  /// @param z New z-axis of the placement.
  DETRAY_HOST
  virtual void add_volume_placement(const dpoint3D<algebra_type> &t,
                                    const dvector3D<algebra_type> &x,
                                    const dvector3D<algebra_type> &z) = 0;

  /// @brief Add surfaces to the volume
  /// @param sf_factory Factory that supplies surfaces for the volume.
  /// @param ctx Geometry context used while adding surfaces.
  DETRAY_HOST
  virtual void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) = 0;

 protected:
  /// Access to builder data
  /// @returns The surface lookup container.
  /// @{
  virtual typename detector_t::surface_lookup_container &surfaces() = 0;
  /// Returns mutable access to the builder's transforms.
  /// @returns The transform container.
  virtual typename detector_t::transform_container &transforms() = 0;
  /// Returns mutable access to the builder's masks.
  /// @returns The mask container.
  virtual typename detector_t::mask_container &masks() = 0;
  /// Returns the hash maps of the detector data that the builder uses for
  /// the deduplication: Either the global cache of the detector builder or
  /// the builder's own cache (standalone builder).
  /// @returns The deduplication cache.
  virtual detail::deduplication_cache &deduplication_cache() = 0;
  /// @}

  /// Central data insertion and deduplication service
  ///
  /// All builders and decorators use these functions to add their
  /// volume-local data to the detector data collections. This keeps the
  /// deduplication in one place, while the builders do not depend on each
  /// other and can be combined freely.
  ///
  /// The data is only deduplicated, if a deduplication cache is given (see
  /// @c deduplication_cache()). The detector builder gives one global cache to all
  /// volume builders on construction. A standalone volume builder owns its
  /// cache.
  ///
  /// Every builder provides the hash function for its own data. The
  /// deduplication cache uses the hash to find candidates and the equality
  /// operator of the data type to confirm them. The hash function has to be
  /// consistent with the equality operator: Identical objects must have the
  /// same hash value.
  /// @{

  /// Add an object to a detector data collection
  ///
  /// @note The collection is accessed with the default geometry context,
  /// since the detector is always built in this context.
  ///
  /// @tparam coll_t the type of the detector collection
  /// @tparam obj_t the type of the object
  /// @tparam hasher_t the type of the hash function for the objects
  ///
  /// @param coll the detector collection that receives the object
  /// @param obj the object that should be added
  /// @param hasher the hash function for the objects
  /// @param opt the builder options
  ///
  /// @returns the index of the object in the detector collection
  template <typename coll_t, typename obj_t, typename hasher_t>
  DETRAY_HOST auto add_to_detector_collection(coll_t &coll, const obj_t &obj,
                                              const hasher_t &hasher,
                                              const volume_builder_options &opt)
      -> dindex {
    if (opt.deduplicate()) {
      if (const dindex idx{this->deduplication_cache().find(coll, obj, hasher)};
          idx != dindex_invalid) {
        DETRAY_DEBUG_HOST("Found identical object at index "
                          << idx << ". Deduplicating...");
        return idx;
      }
    }

    const auto n_objs{static_cast<dindex>(coll.size())};
    coll.push_back(obj);

    return n_objs;
  }

  /// Add all objects of a volume-local collection to a detector collection
  ///
  /// @note The collections are accessed with the default geometry context,
  /// since the detector is always built in this context.
  ///
  /// @tparam coll_t the type of the detector collection
  /// @tparam src_coll_t the type of the volume-local collection
  /// @tparam hasher_t the type of the hash function for the objects
  ///
  /// @param coll the detector collection that receives the objects
  /// @param src the volume-local objects
  /// @param hasher the hash function for the objects
  /// @param opt the builder options
  ///
  /// @returns the detector collection index for every local object
  template <typename coll_t, concepts::range obj_t, typename hasher_t>
  DETRAY_HOST auto add_to_detector_collection(coll_t &coll,
                                              const std::vector<obj_t> &src,
                                              const hasher_t &hasher,
                                              const volume_builder_options &opt)
      -> std::vector<dindex> {
    const auto n_src{static_cast<dindex>(src.size())};

    std::vector<dindex> global_idx{};
    global_idx.reserve(n_src);
    for (dindex i = 0u; i < n_src; ++i) {
      global_idx.push_back(
          add_to_detector_collection(coll, src.at(i), hasher, opt));
    }

    return global_idx;
  }

  /// Add the block of objects that a link refers to from a volume-local
  /// collection to a detector collection
  ///
  /// A block is only reused, if all of its objects are found in the same
  /// order in the detector collection.
  ///
  /// @tparam coll_t the type of the collections
  /// @tparam link_t the type of the link (single index or index range)
  /// @tparam hasher_t the type of the hash function for the objects
  ///
  /// @param coll the detector collection that receives the objects
  /// @param src the volume-local objects
  /// @param link the link into the volume-local collection
  /// @param hasher the hash function for the objects
  /// @param opt the builder options
  ///
  /// @returns the link into the detector collection
  template <typename coll_t, typename link_t, typename hasher_t>
  DETRAY_HOST auto add_linked_objects(coll_t &coll, const coll_t &src,
                                      link_t link, const hasher_t &hasher,
                                      const volume_builder_options &opt)
      -> link_t {
    auto local_idx{link.index()};
    using index_t = decltype(local_idx);

    dindex first{0u};
    dindex n_block{1u};
    if constexpr (concepts::interval<index_t>) {
      first = static_cast<dindex>(local_idx.lower());
      n_block = static_cast<dindex>(local_idx.size());
    } else {
      first = static_cast<dindex>(local_idx);
    }

    const auto n_objs{static_cast<dindex>(coll.size())};
    dindex new_first{n_objs};
    if (opt.deduplicate()) {
      // Candidates for the block start: Objects identical to the first one
      for (const dindex i :
           this->deduplication_cache().find_all(coll, src.at(first), hasher)) {
        if (i + n_block > n_objs) {
          break;
        }
        bool is_match{true};
        for (dindex j = 1u; j < n_block; ++j) {
          if (!(src.at(first + j) == coll.at(i + j))) {
            is_match = false;
            break;
          }
        }
        if (is_match) {
          DETRAY_DEBUG_HOST("Found identical object(s) at index "
                            << i << ". Deduplicating...");
          new_first = i;
          break;
        }
      }
    }

    // No identical block was found: Append the objects
    if (new_first == n_objs) {
      for (dindex j = 0u; j < n_block; ++j) {
        coll.push_back(src.at(first + j));
      }
    }

    // Move the link to the new position
    if constexpr (concepts::interval<index_t>) {
      local_idx.shift(static_cast<long>(new_first) - static_cast<long>(first));
      link.set_index(local_idx);
    } else {
      link.set_index(static_cast<index_t>(new_first));
    }

    return link;
  }
  /// @}
};

/// @brief Decorator for the volume builder.
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
  /// @param ctx Geometry context used while building.
  /// @returns Pointer to the built volume.
  DETRAY_HOST
  auto build(detector_t &det, const volume_builder_options &opt,
             typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type * override {
    return m_builder->build(det, opt, ctx);
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
  /// Access to underlying builder
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
  /// Returns the wrapped builder's deduplication cache.
  /// @returns The deduplication cache.
  detail::deduplication_cache &deduplication_cache() override {
    return m_builder->deduplication_cache();
  }
  /// @}

 private:
  std::unique_ptr<volume_builder_interface<detector_t>> m_builder;
};

}  // namespace detray
