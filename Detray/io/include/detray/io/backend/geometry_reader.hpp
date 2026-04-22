// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <map>
#include <string_view>
#include <utility>
#include <vector>

namespace std {

// Specialize std::hash directly to the shape id type
template <>
struct hash<detray::io::shape_id> {
  auto operator()(detray::io::shape_id m) const noexcept {
    return std::hash<std::size_t>()(static_cast<std::size_t>(m));
  }
};

}  // namespace std

namespace detray::io {

/// @brief Tracking geometry reader backend
///
/// Fills a @c detector_builder from a @c detector_geometry_payload
class geometry_reader {
  /// IO shape ids do not need to coincide with the detector mask ids,
  /// they are shared with ACTS
  using io_shape_id = io::shape_id;

 public:
  /// Tag the reader as "geometry"
  static constexpr std::string_view tag = "geometry";

  /// Payload type that the reader processes
  using payload_type = detector_payload;

  /// Convert a detector @param det from its io payload @param det_data
  template <class detector_t>
  static void from_payload(detector_builder<typename detector_t::metadata,
                                            volume_builder>& det_builder,
                           const payload_type& det_data) {
    DETRAY_VERBOSE_HOST("Reading payload object...");

    // Can hold all types of surface fatory needed for the detector
    using sf_factory_ptr_t =
        std::shared_ptr<surface_factory_interface<detector_t>>;

    // Map the io::shape_id in the payload to the shape types known by
    // the detector under construction
    using shape_registry_t = io::detail::filtered_shape_registry<detector_t>;
    static_assert(
        (!types::contains<shape_registry_t, io::detail::unknown_type> &&
         types::size<shape_registry_t> > 0) ||
        (types::contains<shape_registry_t, io::detail::unknown_type> &&
         types::size<shape_registry_t> > 1));

    DETRAY_DEBUG_HOST("Have " << det_data.volumes.size() << " input volumes");

    // Convert the volumes one-by-one
    for (const auto& vol_data : det_data.volumes) {
      // Get a generic volume builder first and decorate it later
      DETRAY_DEBUG_HOST("Configuring detector builder with new volume '"
                        << vol_data.name << "' of type " << vol_data.type);
      auto vbuilder = det_builder.new_volume(vol_data.type);

      // Set the volume name
      vbuilder->set_name(vol_data.name);

      DETRAY_DEBUG_HOST("Volume placement for " << vol_data.name << " is\n"
                                                << vol_data.transform);
      // Volume placement
      vbuilder->add_volume_placement(
          from_payload<detector_t>(vol_data.transform));

      // Prepare the surface factories (one per shape and surface type)
      std::map<io_shape_id, sf_factory_ptr_t> pt_factories;
      std::map<io_shape_id, sf_factory_ptr_t> sf_factories;

      // Add the surfaces to the factories
      DETRAY_DEBUG_HOST("Discovered " << vol_data.surfaces.size()
                                      << " surfaces in volume payload");
      for (const auto& [sf_idx, sf_data] :
           detray::views::enumerate(vol_data.surfaces)) {
        const std::vector<mask_payload>& mask_data = sf_data.masks;

        // Get a reference to the correct factory collection
        auto& factories{sf_data.type == surface_id::e_portal ? pt_factories
                                                             : sf_factories};

        // @TODO use portal cylinders until intersectors are fixed
        auto shape_id{mask_data.front().shape == io_shape_id::cylinder2
                          ? io_shape_id::portal_cylinder2
                          : mask_data.front().shape};

        if (shape_id == io_shape_id::unknown) {
          std::string err_str{"Unknown mask shape id in geometry file!"};
          DETRAY_FATAL_HOST(err_str);
          throw std::invalid_argument(err_str);
        }

        // Check if a fitting factory already exists. If not, add it
        // dynamically
        const auto key{shape_id};
        if (auto search = factories.find(key); search == factories.end()) {
          DETRAY_DEBUG_HOST(
              "Creating new surface factory for shape id: " << key);
          factories[key] =
              std::move(init_factory<shape_registry_t, detector_t>(shape_id));
        }

        DETRAY_DEBUG_HOST("-> Surface #" << sf_idx << " is " << sf_data.type);
        DETRAY_DEBUG_HOST("--> shape = " << shape_id);
        for ([[maybe_unused]] const auto& m_data : mask_data) {
          DETRAY_DEBUG_HOST("--> boundaries = ["
                            << DETRAY_LOG_VECTOR(m_data.boundaries) << "]");
        }
        DETRAY_DEBUG_HOST("--> pos = ["
                          << DETRAY_LOG_VECTOR(sf_data.transform.tr) << "]");

        // Add the data to the factory
        factories.at(key)->push_back(from_payload<detector_t>(sf_data));
      }

      // Add all portals and surfaces to the volume
      typename detector_t::geometry_context geo_ctx{};
      for (auto [key, pt_factory] : pt_factories) {
        DETRAY_DEBUG_HOST("Adding '" << key
                                     << "' portal factory to volume builder");
        vbuilder->add_surfaces(pt_factory, geo_ctx);
      }

      for (auto [key, sf_factory] : sf_factories) {
        DETRAY_DEBUG_HOST("Adding '" << key
                                     << "' surface factory to volume builder");
        vbuilder->add_surfaces(sf_factory, geo_ctx);
      }
    }

    // @TODO: Implement voume finder IO
  }

  /// @returns a surface transform from its io payload @param trf_data
  template <class detector_t>
  static typename detector_t::transform3_type from_payload(
      const transform_payload& trf_data) {
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using vector3_t = dvector3D<algebra_t>;

    vector3_t t{static_cast<scalar_t>(trf_data.tr[0]),
                static_cast<scalar_t>(trf_data.tr[1]),
                static_cast<scalar_t>(trf_data.tr[2])};
    vector3_t x{static_cast<scalar_t>(trf_data.rot[0]),
                static_cast<scalar_t>(trf_data.rot[1]),
                static_cast<scalar_t>(trf_data.rot[2])};
    vector3_t y{static_cast<scalar_t>(trf_data.rot[3]),
                static_cast<scalar_t>(trf_data.rot[4]),
                static_cast<scalar_t>(trf_data.rot[5])};
    vector3_t z{static_cast<scalar_t>(trf_data.rot[6]),
                static_cast<scalar_t>(trf_data.rot[7]),
                static_cast<scalar_t>(trf_data.rot[8])};

    return dtransform3D<algebra_t>{t, x, y, z};
  }

  /// @returns surface data for a surface factory from a surface io payload
  /// @param trf_data
  template <class detector_t>
  static surface_data<detector_t> from_payload(const surface_payload& sf_data) {
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    using scalar_t = dscalar<typename detector_t::algebra_type>;

    // Transcribe mask boundaries onto correct vector type
    std::vector<std::vector<scalar_t>> mask_boundaries{};
    for (const auto& mask_data : sf_data.masks) {
      std::vector<scalar_t>& boundaries = mask_boundaries.emplace_back();
      std::ranges::copy(mask_data.boundaries, std::back_inserter(boundaries));
    }

    // If the concentric cylinder is shifted in z, discard the shift
    // and put it in the mask boundaries instead. Everything else
    // will be ignored
    // @TODO: Remove this for cylinders again once 2 solution intersection
    // work
    auto trf = from_payload<detector_t>(sf_data.transform);
    if (sf_data.masks.front().shape == io_shape_id::portal_cylinder2 ||
        sf_data.masks.front().shape == io_shape_id::cylinder2) {
      const auto z_shift{static_cast<scalar_t>(trf.translation()[2])};

      for (auto& mask_boundary : mask_boundaries) {
        mask_boundary[concentric_cylinder2D::e_lower_z] += z_shift;
        mask_boundary[concentric_cylinder2D::e_upper_z] += z_shift;
      }

      // Set the transform to identity afterwards
      trf = decltype(trf){};
    }

    const std::size_t sf_idx{
        sf_data.index_in_coll.has_value()
            ? *(sf_data.index_in_coll)
            : detray::detail::invalid_value<std::size_t>()};

    std::vector<nav_link_t> vol_links{};
    vol_links.reserve(sf_data.masks.size());
    for (const auto& mask_pl : sf_data.masks) {
      vol_links.push_back(static_cast<nav_link_t>(
          detail::basic_converter::from_payload(mask_pl.volume_link)));
    }

    assert(mask_boundaries.size() == vol_links.size());

    return {sf_data.type,
            trf,
            std::move(vol_links),
            std::move(mask_boundaries),
            static_cast<dindex>(sf_idx),
            sf_data.source};
  }

 private:
  /// Determines the surface shape type from the IO id @param shape_id in the
  /// payload
  ///
  /// @return the corresponding surface factory.
  template <typename shape_registry_t, typename detector_t, std::size_t I = 0u>
  static std::shared_ptr<surface_factory_interface<detector_t>> init_factory(
      const io_shape_id shape_id) {
    // Get the next mask shape type
    using shape_t = types::at<shape_registry_t, I>;

    // Test whether this shape exists in detector
    if constexpr (!std::is_same_v<shape_t, io::detail::unknown_type>) {
      // Shape index of surface data found
      if (shape_id == types::id<io::shape_registry, shape_t>) {
        return std::make_shared<surface_factory<detector_t, shape_t>>();
      }
    }
    // Test next shape id
    if constexpr (I < types::size<shape_registry_t> - 1u) {
      return init_factory<shape_registry_t, detector_t, I + 1>(shape_id);
    } else {
      std::string err_str{
          "Given shape id could not be matched to a mask type: " +
          std::to_string(static_cast<std::int64_t>(shape_id))};

      DETRAY_FATAL_HOST(err_str);
      throw std::invalid_argument(err_str);
    }

    // Cannot be reached
    return {nullptr};
  }
};

}  // namespace detray::io
