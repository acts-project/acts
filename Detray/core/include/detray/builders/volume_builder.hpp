// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder_interface.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace detray {

namespace detail {

/// A functor to update the mask index in surface descriptors
struct mask_index_update;

}  // namespace detail

/// @brief Provides basic functionality to build detector volumes
template <typename detector_t>
class volume_builder : public volume_builder_interface<detector_t> {
  static_assert(concepts::detector<detector_t>);

 public:
  /// Scalar type used by this detector's geometry.
  using scalar_t = dscalar<typename detector_t::algebra_type>;
  /// Volume descriptor type produced by the builder.
  using volume_type = typename detector_t::volume_type;
  /// Identifier enum for the detector's geometry object categories.
  using geo_obj_ids = typename detector_t::geo_obj_ids;

  /// Parametrized Constructor
  ///
  /// @param id flags the type of volume geometry (e.g. cylindrical, cuboid)
  /// @param idx the index of the volume in the detector volume container
  explicit volume_builder(const volume_id id, const dindex idx = 0)
      : m_volume{id} {
    m_volume.set_index(idx);
    m_volume.set_material(volume_type::material_id::e_none, dindex_invalid);

    // The first acceleration data structure in every volume is a brute
    // force method that will at least contain the portals
    m_volume.template set_accel_link<
        static_cast<typename volume_type::object_id>(0)>(
        detector_t::accel::id::e_surface_default, 0);

    DETRAY_VERBOSE_HOST("Created builder for volume: " << idx);
  };

  /// @returns the volume index in the detector volume container
  DETRAY_HOST
  auto vol_index() const -> dindex override { return m_volume.index(); }

  /// Toggles whether sensitive surfaces are added to the brute force method
  /// @param toggle True to retain sensitive surfaces in the brute-force accelerator.
  DETRAY_HOST
  void has_accel(bool toggle) override { m_has_accel = toggle; }

  /// @returns whether sensitive surfaces are added to the brute force method
  DETRAY_HOST
  bool has_accel() const override { return m_has_accel; }

  /// Sets the name @p volume_name for the volume
  /// @param volume_name Name to assign to the volume.
  DETRAY_HOST void set_name(std::string volume_name) override {
    m_volume_name = std::move(volume_name);
    DETRAY_VERBOSE_HOST("Set volume name: " << m_volume_name);
  }

  /// @returns the name of the volume
  DETRAY_HOST std::string_view name() override {
    if (m_volume_name.empty()) {
      // Consistent default after volume index is known
      m_volume_name = "unknown(volume_" + std::to_string(vol_index()) + ")";
    }
    return m_volume_name;
  }

  /// Access to the volume under construction - const
  /// @returns Const reference to the volume descriptor.
  DETRAY_HOST
  auto operator()() const -> const typename detector_t::volume_type& override {
    return m_volume;
  }

  /// Access to the volume under construction - non-const
  /// @returns Mutable reference to the volume descriptor.
  DETRAY_HOST
  auto operator()() -> typename detector_t::volume_type& override {
    return m_volume;
  }

  /// Build the volume with internal surfaces and portals and add it to the
  /// detector instance @p det
  /// @param det Detector that receives the built geometry.
  /// @param opt Options controlling volume construction.
  /// @param ctx Geometry context used when storing transforms.
  /// @returns Pointer to the volume descriptor added to the detector.
  DETRAY_HOST
  auto build(detector_t& det, const volume_builder_options& opt,
             typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type* override {
    DETRAY_VERBOSE_HOST("Build surfaces...");

    assert(!m_surfaces.empty());
    assert(!m_transforms.empty());
    assert(!m_masks.all_empty());

    // Prepare volume data
    m_volume.set_index(static_cast<dindex>(det.volumes().size()));

    m_volume.set_transform(add_transform(det, m_trf, opt.deduplicate(), ctx));

    // Add all data from the builder to the detector containers
    add_to_detector(opt, ctx, det);

    // Reset after the data was added to the detector
    m_surfaces.clear();
    m_transforms.clear(ctx);
    m_masks.clear_all();

    DETRAY_VERBOSE_HOST("Successfully built "
                        << m_volume.n_surfaces()
                        << " surfaces for volume: " << this->name());

    // Pass to decorator builders
    return &(det._volumes.back());
  }

  /// Adds a placement transform @param trf for the volume
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::transform3_type& trf = {}) override {
    m_trf = trf;
  }

  /// Constructs a placement transform with identity rotation and translation
  /// @param t for the volume
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type& t) override {
    m_trf = typename detector_t::transform3_type{t};
  }

  /// Constructs a placement transform from axes @p x and @p z
  /// and the translation @p t for the volume
  /// @param t Translation of the volume in global coordinates.
  /// @param x Local x axis in global coordinates.
  /// @param z Local z axis in global coordinates.
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type& t,
      const typename detector_t::vector3_type& x,
      const typename detector_t::vector3_type& z) override {
    m_trf = typename detector_t::transform3_type{t, z, x};
  }

  /// Add data for (a) new surface(s) to the builder
  /// @param sf_factory Factory that creates the surfaces and masks.
  /// @param ctx Geometry context used when constructing transforms.
  DETRAY_HOST
  void add_surfaces(
      std::shared_ptr<surface_factory_interface<detector_t>> sf_factory,
      typename detector_t::geometry_context ctx = {}) override {
    DETRAY_VERBOSE_HOST("Add surface factory:");

    (*sf_factory)(m_volume, m_surfaces, m_transforms, m_masks, ctx);
  }

 protected:
  /// @returns Access to the surface descriptor data
  typename detector_t::surface_lookup_container& surfaces() override {
    return m_surfaces;
  }

  /// @returns Access to the surface/volume transform data
  typename detector_t::transform_container& transforms() override {
    return m_transforms;
  }

  /// @returns Access to the surface mask data
  typename detector_t::mask_container& masks() override { return m_masks; }

  /// Adds a new full set of volume components (e.g. transforms or masks)
  /// to the global detector data stores and updates all links.
  ///
  /// @param opt Options controlling deduplication and geometry storage.
  /// @param ctx is the geometry_context of the call
  /// @param det is the detector instance that the volume should be added to
  ///
  /// @note can throw an exception if input data is inconsistent
  template <geo_obj_ids surface_id = static_cast<geo_obj_ids>(0)>
  DETRAY_HOST auto add_to_detector(
      const volume_builder_options& opt,
      const typename detector_t::geometry_context ctx,
      detector_t& det) noexcept(false) -> void {
    const bool deduplicate{opt.deduplicate()};

    // Global transform index for every volume local transform
    std::vector<dindex> global_trf_idx{};
    // Whether a global transform is used by more than one surface
    std::vector<bool> is_shared_trf{};

    const auto trf_offset = det.transform_store().size(ctx);
    if (deduplicate) {
      DETRAY_VERBOSE_HOST("-> Deduplicate transforms and masks");

      // Check every local transform against the detector transforms
      // (Any duplication within the new data is resolved the same way)
      const auto n_local_trfs{m_transforms.size(ctx)};
      global_trf_idx.reserve(n_local_trfs);
      std::vector<dindex> found_trfs{};
      for (dindex i = 0u; i < n_local_trfs; ++i) {
        const auto n_trfs{det.transform_store().size(ctx)};
        const dindex trf_idx{
            add_transform(det, m_transforms.at(i, ctx), true, ctx)};
        // An identical transform was already in the detector store
        if (trf_idx < n_trfs) {
          found_trfs.push_back(trf_idx);
        }
        global_trf_idx.push_back(trf_idx);
      }

      // Flag the transforms that are used by more than one local entry or
      // that were in the detector before
      is_shared_trf.resize(det.transform_store().size(ctx), false);
      for (const dindex trf_idx : found_trfs) {
        is_shared_trf.at(trf_idx) = true;
      }

      // Add the masks
      add_masks(det, std::make_index_sequence<detector_t::masks::n_types>{});
    } else {
      // Append all transforms
      det._transforms.append(std::move(m_transforms), ctx);
    }

    // Surface index offset in the global detector container
    auto sf_offset{static_cast<dindex>(det.surfaces().size())};

    /// Find the surface range specified by @param sf_id
    auto find_range = [&](auto sf_id) {
      // Compare given id to surface identifier
      auto is_sf_type = [sf_id](const auto& sf) { return sf.id() == sf_id; };

      auto first = static_cast<dindex>(
          math::abs(std::ranges::find_if(m_surfaces, is_sf_type) -
                    std::begin(m_surfaces)));

      auto last = static_cast<dindex>(math::abs(
          std::ranges::rend(m_surfaces) -
          std::ranges::find_if(std::ranges::rbegin(m_surfaces),
                               std::ranges::rend(m_surfaces), is_sf_type)));

      // Set correct empty range, otherwise shift by global surface offset
      return (first >= last)
                 ? dindex_range{}
                 : dindex_range{first + sf_offset, last + sf_offset};
    };

    m_volume.template update_sf_link<surface_id::e_portal>(
        find_range(surface_id::e_portal));

    m_volume.template update_sf_link<surface_id::e_sensitive>(
        find_range(surface_id::e_sensitive));

    m_volume.template update_sf_link<surface_id::e_passive>(
        find_range(surface_id::e_passive));

    // Update mask and transform index of surfaces and set the
    // correct index of the surface in container
    std::size_t n_portals{0u};
    // Count the surfaces that use a given transform (only if deduplicated)
    std::vector<dindex> n_sf_per_trf(is_shared_trf.size(), 0u);
    for (auto& sf_desc : m_surfaces) {
      if (deduplicate) {
        // Mask links were already updated during the deduplication
        const dindex trf_idx{global_trf_idx.at(sf_desc.transform())};
        sf_desc.set_transform(trf_idx);
        ++n_sf_per_trf.at(trf_idx);
      } else {
        det._masks.template visit<detail::mask_index_update>(sf_desc.mask(),
                                                             sf_desc);
        sf_desc.update_transform(trf_offset);
      }
      sf_desc.set_volume(m_volume.index());
      sf_desc.set_index(sf_offset++);

      // Sensitive surfaces with a shared transform cannot be aligned
      // independently
      if (deduplicate && sf_desc.is_sensitive() &&
          (is_shared_trf.at(sf_desc.transform()) ||
           n_sf_per_trf.at(sf_desc.transform()) > 1u)) {
        warn_shared_sensitive_trf(det, sf_desc);
      }

      if (sf_desc.is_portal()) {
        ++n_portals;
      }

      det._surfaces.insert(sf_desc);
    }

    // Place the appropriate surfaces in the brute force search method.
    constexpr auto default_acc_id{detector_t::accel::id::e_surface_default};

    // Strip the source link from the lookup data structure
    typename detector_t::surface_container descriptors;
    descriptors.reserve(m_surfaces.size());
    std::ranges::transform(
        m_surfaces, std::back_inserter(descriptors),
        [](typename detector_t::surface_lookup_container::value_type& sf) {
          return static_cast<
              typename detector_t::surface_container::value_type>(sf);
        });

    // Add portals to brute force navigation method
    if (m_has_accel) {
      DETRAY_VERBOSE_HOST("-> Volume has acceleration structure:");

      typename detector_t::surface_container portals{};
      portals.reserve(n_portals);

      std::ranges::copy_if(
          descriptors, std::back_inserter(portals),
          [](auto& sf_desc) { return !sf_desc.is_sensitive(); });

      // Add only the portals to the brute force method
      DETRAY_VERBOSE_HOST(
          "-> Register only portals/passives with brute force "
          "accelerator");

      det._accelerators.template push_back<default_acc_id>(std::move(portals));
    } else {
      DETRAY_VERBOSE_HOST(
          "-> Register all surfaces with brute force accelerator");

      // Add all surfaces to the brute force method
      det._accelerators.template push_back<default_acc_id>(
          std::move(descriptors));
    }

    m_volume.template set_accel_link<surface_id>(
        default_acc_id,
        det.accelerator_store().template size<default_acc_id>() - 1u);

    // Append masks
    if (!deduplicate) {
      det._masks.append(std::move(m_masks));
    }

    // Finally, add the volume descriptor to the detector
    det._volumes.push_back(m_volume);
  }

 private:
  /// Add a transform to the detector transform store
  ///
  /// @param det the detector that receives the transform
  /// @param trf the transform that should be added
  /// @param deduplicate whether to reuse an identical transform in @p det
  /// @param ctx the geometry context of the call
  ///
  /// @returns the index of the transform in the detector transform store
  DETRAY_HOST auto add_transform(
      detector_t& det, const typename detector_t::transform3_type& trf,
      const bool deduplicate,
      const typename detector_t::geometry_context ctx) const -> dindex {
    const auto n_trfs{static_cast<dindex>(det.transform_store().size(ctx))};

    // Test only against the transforms that are already in the detector
    if (deduplicate) {
      for (dindex i = 0u; i < n_trfs; ++i) {
        if (trf == det.transform_store().at(i, ctx)) {
          DETRAY_DEBUG_HOST("Found identical transform at index "
                            << i << ". Deduplicating...");
          return i;
        }
      }
    }

    det._transforms.push_back(trf, ctx);

    return n_trfs;
  }

  /// Add the masks of all mask types to the detector and deduplicate them
  ///
  /// @tparam I the indices of the mask types in the mask store
  ///
  /// @param det the detector that receives the masks
  template <std::size_t... I>
  DETRAY_HOST void add_masks(detector_t& det, std::index_sequence<I...>) {
    (add_masks<static_cast<typename detector_t::masks::id>(I)>(det), ...);
  }

  /// Add the masks of type @tparam mask_id to the detector and reuse
  /// identical masks that are already in the detector store
  ///
  /// The masks of a surface are compared as a block, i.e. a mask range is
  /// only reused if all masks of the surface are found in the same
  /// order in the detector store. The surface mask links are updated.
  ///
  /// @param det the detector that receives the masks
  template <typename detector_t::masks::id mask_id>
  DETRAY_HOST void add_masks(detector_t& det) {
    const auto& local_coll = m_masks.template get<mask_id>();
    if (local_coll.empty()) {
      return;
    }

    // non-const access
    auto& det_coll = det._masks.template get<mask_id>();
    [[maybe_unused]] const auto n_masks_before{det_coll.size()};

    for (auto& sf_desc : m_surfaces) {
      auto mask_link{sf_desc.mask()};
      if (mask_link.id() != mask_id) {
        continue;
      }

      // The volume local masks of the surface
      auto local_idx{mask_link.index()};
      dindex first{0u};
      dindex n_sf_masks{1u};
      if constexpr (concepts::interval<decltype(local_idx)>) {
        first = static_cast<dindex>(local_idx.lower());
        n_sf_masks = static_cast<dindex>(local_idx.size());
      } else {
        first = static_cast<dindex>(local_idx);
      }

      // Search the detector for the same sequence of masks
      const auto n_det_masks{static_cast<dindex>(det_coll.size())};
      dindex new_first{n_det_masks};
      for (dindex i = 0u; i + n_sf_masks <= n_det_masks; ++i) {
        bool is_match{true};
        for (dindex j = 0u; j < n_sf_masks; ++j) {
          if (!(local_coll.at(first + j) == det_coll.at(i + j))) {
            is_match = false;
            break;
          }
        }
        if (is_match) {
          DETRAY_DEBUG_HOST("Found identical mask(s) at index "
                            << i << ". Deduplicating...");
          new_first = i;
          break;
        }
      }

      // No duplicate was found, append the new mask(s)
      if (new_first == n_det_masks) {
        for (dindex j = 0u; j < n_sf_masks; ++j) {
          det_coll.push_back(local_coll.at(first + j));
        }
      }

      // Update the surface mask link
      if constexpr (concepts::interval<decltype(local_idx)>) {
        local_idx.shift(static_cast<long>(new_first) -
                        static_cast<long>(first));
        mask_link.set_index(local_idx);
      } else {
        mask_link.set_index(new_first);
      }
      sf_desc.set_mask(mask_link);
    }

    DETRAY_VERBOSE_HOST("-> Appended " << det_coll.size() - n_masks_before
                                       << " of " << local_coll.size()
                                       << " masks of type " << mask_id
                                       << " to detector masks");
  }

  /// Warn if the sensitive surface @param sf_desc shares its transform
  /// with another sensitive surface in the detector @param det
  ///
  /// @note Sensitive surfaces with a shared transform cannot be misaligned
  /// independently
  DETRAY_HOST void warn_shared_sensitive_trf(
      const detector_t& det,
      const typename detector_t::surface_lookup_container::value_type& sf_desc)
      const {
    for (const auto& other : det.surfaces()) {
      if (other.is_sensitive() && other.transform() == sf_desc.transform()) {
        DETRAY_WARN_HOST(m_volume_name
                         << ": Sensitive surfaces share the transform "
                         << sf_desc.transform() << " (surfaces " << other.index()
                         << " and " << sf_desc.index()
                         << "). They cannot be aligned independently.");
        return;
      }
    }
  }

  /// Whether the volume will get an acceleration structure
  bool m_has_accel{false};

  /// The name of the volume
  std::string m_volume_name{};

  /// Volume descriptor of the volume under construction
  typename detector_t::volume_type m_volume{};
  /// Placement of the volume under construction
  typename detector_t::transform3_type m_trf{};

  /// Data of contained surfaces
  /// @{
  typename detector_t::surface_lookup_container m_surfaces{};
  typename detector_t::transform_container m_transforms{};
  typename detector_t::mask_container m_masks{};
  /// @}
};

namespace detail {

/// A functor to update the mask index in surface objects
struct mask_index_update {
  /// Update a surface's mask index after mask groups are appended.
  /// @param group Mask group whose size determines the new offset.
  /// @param sf Surface descriptor whose mask link is updated.
  template <typename group_t, typename index_t, typename surface_t>
  DETRAY_HOST inline void operator()(const group_t& group,
                                     const index_t& /*index*/,
                                     surface_t& sf) const {
    sf.update_mask(static_cast<dindex>(std::size(group)));
  }
};

}  // namespace detail

}  // namespace detray
