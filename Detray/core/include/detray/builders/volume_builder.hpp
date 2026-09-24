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
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <algorithm>
#include <memory>
#include <sstream>
#include <string>

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
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using surface_type = typename detector_t::surface_type;
  using volume_type = typename detector_t::volume_type;
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
  DETRAY_HOST
  void has_accel(bool toggle) override { m_has_accel = toggle; }

  /// @returns whether sensitive surfaces are added to the brute force method
  DETRAY_HOST
  bool has_accel() const override { return m_has_accel; }

  /// Sets the name @param volume_name for the volume
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
  DETRAY_HOST
  auto operator()() const -> const typename detector_t::volume_type& override {
    return m_volume;
  }

  /// Access to the volume under construction - non-const
  DETRAY_HOST
  auto operator()() -> typename detector_t::volume_type& override {
    return m_volume;
  }

  /// Build the volume with internal surfaces and portals and add it to the
  /// detector instance @param det
  DETRAY_HOST
  auto build(detector_t& det, typename detector_t::geometry_context ctx = {}) ->
      typename detector_t::volume_type* override {
    DETRAY_VERBOSE_HOST("Build surfaces...");

    assert(!m_surfaces.empty());
    assert(!m_transforms.empty());
    assert(!m_masks.all_empty());

    // Prepare volume data
    m_volume.set_index(static_cast<dindex>(det.volumes().size()));

    m_volume.set_transform(det.transform_store().size());
    det._transforms.push_back(m_trf);

    // Add all data from the builder to the detector containers
    add_to_detector(ctx, det);

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

  /// Constructs a placement transform from axes @param x and @param z
  /// and the translation @param t for the volume
  DETRAY_HOST
  void add_volume_placement(
      const typename detector_t::point3_type& t,
      const typename detector_t::vector3_type& x,
      const typename detector_t::vector3_type& z) override {
    m_trf = typename detector_t::transform3_type{t, z, x};
  }

  /// Add data for (a) new surface(s) to the builder
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
  /// @param ctx is the geometry_context of the call
  /// @param det is the detector instance that the volume should be added to
  ///
  /// @note can throw an exception if input data is inconsistent
  template <geo_obj_ids surface_id = static_cast<geo_obj_ids>(0)>
  DETRAY_HOST auto add_to_detector(
      const typename detector_t::geometry_context ctx,
      detector_t& det) noexcept(false) -> void {
    // Append transforms
    const auto trf_offset = det.transform_store().size(ctx);
    det._transforms.append(std::move(m_transforms), ctx);

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

    // Make sure, the portals fit the volume boundaries, otherwise clip them
    switch (m_volume.id()) {
      case volume_id::e_cylinder: {
        using masks = typename detector_t::masks;
        using cylinder_t = mask<concentric_cylinder2D, algebra_t>;
        using disc_t = mask<ring2D, algebra_t>;

        if constexpr (detray::types::contains<masks, cylinder_t> &&
                      detray::types::contains<masks, disc_t>) {
          DETRAY_VERBOSE_HOST("Clipping portals to cylinder volume shape...");

          auto& cyls =
              m_masks.template get<masks::id::e_concentric_cylinder2D>();
          auto& discs = m_masks.template get<masks::id::e_ring2D>();

          constexpr auto inv{detail::invalid_value<scalar_t>()};

          // Find z extent
          scalar_t min_z{inv};
          scalar_t max_z{-inv};
          for (const surface_type& sf_desc : m_surfaces) {
            if (sf_desc.is_portal() &&
                sf_desc.mask().id() == masks::id::e_ring2D) {
              const auto& t =
                  m_transforms.at(sf_desc.transform()).translation();
              min_z = math::min(min_z, t[2]);
              max_z = math::max(max_z, t[2]);
            }
          }
          if (min_z >= max_z) {
            std::stringstream err{};
            err << "Detected invalid cylinder volume extent: min z: " << min_z
                << "mm, max z: " << max_z << "mm";
            DETRAY_FATAL_HOST(err.str());
            throw std::invalid_argument(err.str());
          }

          // Find radial extent and clip cylinders in z
          scalar_t min_r{inv};
          scalar_t max_r{-inv};
          for (cylinder_t& c : cyls) {
            const scalar_t r{c[concentric_cylinder2D::e_r]};
            min_r = math::min(min_r, r);
            max_r = math::max(max_r, r);

            if ((c[concentric_cylinder2D::e_lower_z] < min_z &&
                 c[concentric_cylinder2D::e_upper_z] < min_z) ||
                (c[concentric_cylinder2D::e_lower_z] > max_z &&
                 c[concentric_cylinder2D::e_upper_z] > max_z)) {
              DETRAY_ERROR_HOST("Portal ["
                                << c
                                << "] lies completely outside cylinder volume '"
                                << m_volume_name << "' z: [" << min_z << ", "
                                << max_z << "] and needs to be removed!");
              continue;
            }

            DETRAY_DEBUG_HOST("Cylinder: " << c);
            c[concentric_cylinder2D::e_lower_z] =
                math::max(min_z, c[concentric_cylinder2D::e_lower_z]);
            c[concentric_cylinder2D::e_upper_z] =
                math::min(max_z, c[concentric_cylinder2D::e_upper_z]);
            DETRAY_DEBUG_HOST("-> clipped: " << c);
          }

          // Beampipe or world volume (no inner cylinder, r is exactly eq.)
          if (min_r == max_r) {
            min_r = 0.f;
          }
          if (min_r >= max_r) {
            std::stringstream err{};
            err << "Detected invalid cylinder volume extent: min r: " << min_r
                << "mm, max r: " << max_r << "mm";
            DETRAY_FATAL_HOST(err.str());
            throw std::invalid_argument(err.str());
          }

          // Clip disc portals to cylinder radius
          for (disc_t& d : discs) {
            if ((d[ring2D::e_inner_r] < min_r &&
                 d[ring2D::e_outer_r] < min_r) ||
                (d[ring2D::e_inner_r] > max_r &&
                 d[ring2D::e_outer_r] > max_r)) {
              DETRAY_ERROR_HOST(
                  "Portal ["
                  << d << "] lies completely outside cylinder volume '"
                  << m_volume_name << "' radius: [" << min_r << ", " << max_r
                  << "] and needs to be removed!");
              continue;
            }
            DETRAY_DEBUG_HOST("Disc: " << d);
            d[ring2D::e_inner_r] = math::max(min_r, d[ring2D::e_inner_r]);
            d[ring2D::e_outer_r] = math::min(max_r, d[ring2D::e_outer_r]);
            DETRAY_DEBUG_HOST("-> clipped: " << d);
          }
        } else {
          const std::string err{
              "Detector with cylinder volumes does not contain cylinder and "
              "disc types in metadata!"};
          DETRAY_FATAL_HOST(err);
          throw std::invalid_argument(err);
        }
        break;
      }
      case volume_id::e_rectangle: {
        DETRAY_DEBUG_HOST(
            "Portal clipping not implemented for rectangle volumes");
        break;
      }
      case volume_id::e_trapezoid: {
        DETRAY_DEBUG_HOST(
            "Portal clipping not implemented for trapezoid volumes");
        break;
      }
      case volume_id::e_cone: {
        DETRAY_DEBUG_HOST("Portal clipping not implemented for cone volumes");
        break;
      }
      case volume_id::e_cuboid: {
        DETRAY_DEBUG_HOST("Portal clipping not implemented for cuboid volumes");
        break;
      }
      case volume_id::e_unknown: {
        DETRAY_WARN_HOST("Unknown volume shape: portal clipping impossible");
        break;
      }
      default: {
        const std::string err{"Unknown error during volume portal clipping"};
        DETRAY_FATAL_HOST(err);
        throw std::invalid_argument(err);
      }
    }

    // Update mask and transform index of surfaces and set the
    // correct index of the surface in container
    std::size_t n_portals{0u};
    for (auto& sf_desc : m_surfaces) {
      det._masks.template visit<detail::mask_index_update>(sf_desc.mask(),
                                                           sf_desc);
      sf_desc.set_volume(m_volume.index());
      sf_desc.update_transform(trf_offset);
      sf_desc.set_index(sf_offset++);

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
    det._masks.append(std::move(m_masks));

    // Finally, add the volume descriptor to the detector
    det._volumes.push_back(m_volume);
  }

 private:
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
  template <typename group_t, typename index_t, typename surface_t>
  DETRAY_HOST inline void operator()(const group_t& group,
                                     const index_t& /*index*/,
                                     surface_t& sf) const {
    sf.update_mask(static_cast<dindex>(std::size(group)));
  }
};

}  // namespace detail

}  // namespace detray
