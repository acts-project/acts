// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/grid/concepts.hpp"

// System include(s)
#include <algorithm>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Tracking geometry writer backend
///
/// Fills a @c detector_payload from a @c detector instance
class geometry_writer {
 public:
  /// Tag the writer as "geometry"
  static constexpr std::string_view tag = "geometry";

  /// Payload type that the reader processes
  using payload_type = detector_payload;

  /// Convert the header information into its payload
  template <typename detector_t>
  static geo_header_payload header_to_payload(const detector_t& det,
                                              const std::string_view det_name) {
    geo_header_payload header_data;

    header_data.common = detail::basic_converter::to_payload(det_name, tag);

    header_data.sub_header.emplace();
    auto& geo_sub_header = header_data.sub_header.value();
    geo_sub_header.n_volumes = det.volumes().size();
    geo_sub_header.n_surfaces = det.surfaces().size();

    return header_data;
  }

  /// Convert a detector @param det into its io payload
  template <typename detector_t>
  static payload_type to_payload(const detector_t& det,
                                 const typename detector_t::name_map& names) {
    payload_type det_data;
    det_data.volumes.reserve(det.volumes().size());

    for (const auto& vol : det.volumes()) {
      if (names.contains(vol.index())) {
        det_data.volumes.push_back(to_payload(vol, det, names.at(vol.index())));
      } else {
        det_data.volumes.push_back(to_payload(vol, det, ""));
      }
    }

    return det_data;
  }

  /// Convert a surface transform @param trf into its io payload
  template <typename detector_t>
  static transform_payload to_payload(
      const typename detector_t::transform3_type& trf) {
    transform_payload trf_data;

    const auto& t = trf.translation();
    const auto& x = trf.x();
    const auto& y = trf.y();
    const auto& z = trf.z();

    trf_data.tr = {t[0], t[1], t[2]};
    trf_data.rot = {x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]};

    return trf_data;
  }

  /// Convert a surface mask @param m into its io payload
  template <typename mask_t>
    requires(!std::is_same_v<typename mask_t::shape, void>)
  static mask_payload to_payload(const mask_t& m) {
    mask_payload mask_data;

    mask_data.shape = io::detail::get_id<typename mask_t::shape>();

    mask_data.volume_link =
        detail::basic_converter::to_payload(m.volume_link());

    mask_data.boundaries.resize(mask_t::boundaries::e_size);
    std::ranges::copy(m.values(), std::begin(mask_data.boundaries));

    return mask_data;
  }

  /// Convert a detector surface @param sf into its io payload
  template <typename detector_t>
  static surface_payload to_payload(const geometry::surface<detector_t>& sf,
                                    std::size_t sf_idx) {
    using algebra_t = typename detector_t::algebra_type;
    surface_payload sf_data;

    sf_data.index_in_coll = sf_idx;
    sf_data.type = sf.id();
    sf_data.identifier = sf.identifier().value();
    sf_data.transform = to_payload<detector_t>(sf.transform({}));
    sf_data.masks = sf.template visit_mask<get_mask_link_payload>();
    if (sf.has_material()) {
      sf_data.material =
          sf.template visit_material<get_material_link_payload<algebra_t>>();
    }
    sf_data.source = sf.source();

    return sf_data;
  }

  /// Convert a detector portal @param sf into its io payload
  template <typename detector_t>
  static volume_payload to_payload(
      const typename detector_t::volume_type& vol_desc, const detector_t& det,
      const std::string_view name) {
    volume_payload vol_data;

    vol_data.index = detail::basic_converter::to_payload(vol_desc.index());
    vol_data.name = name;
    vol_data.transform =
        to_payload<detector_t>(det.transform_store().at(vol_desc.transform()));
    vol_data.type = vol_desc.id();

    // Count the surfaces belonging to this volume
    std::size_t sf_idx{0};

    for (const auto& sf_desc : det.surfaces()) {
      if (sf_desc.volume() == vol_desc.index()) {
        vol_data.surfaces.push_back(
            to_payload(geometry::surface{det, sf_desc}, sf_idx++));
      }
    }

    // Only run the query, if object type is contained in volume
    const auto& link = vol_desc.accel_link();
    // Initialize the std::optional
    if (link.size() > 1u) {
      vol_data.acc_links.emplace();
    }
    // Skip the first acceleration structure which exists in every volume
    // and is handled automatically during detector building
    for (unsigned int i = 1u; i < link.size(); ++i) {
      const auto& l = link[i];
      if (!l.is_invalid()) {
        const auto aclp =
            det.accelerator_store().template visit<get_acc_link_payload>(l);
        vol_data.acc_links->push_back(aclp);
      }
    }

    return vol_data;
  }

 private:
  /// Retrieve @c mask_payload from mask_store element
  struct get_mask_link_payload {
    template <typename mask_group_t, typename idx_range_t>
    inline auto operator()(const mask_group_t& mask_group,
                           const idx_range_t& idx_range) const {
      std::vector<mask_payload> payloads{};

      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        payloads.push_back(geometry_writer::to_payload(mask));
      }

      return payloads;
    }
  };

  /// Retrieve @c material_link_payload from material_store element
  template <detray::concepts::algebra algebra_t>
  struct get_material_link_payload {
    template <typename material_group_t, typename index_t>
    constexpr auto operator()(const material_group_t& /*unused*/,
                              const index_t& index) const {
      using material_t = typename material_group_t::value_type;

      // Find the correct material type index
      if constexpr (detray::concepts::grid<material_t>) {
        return detail::basic_converter::to_payload(
            io::detail::get_id<material_t>(), index);
      } else {
        return detail::basic_converter::to_payload(
            io::detail::get_id<algebra_t, material_t>(), index);
      }
    }
  };

  /// Retrieve @c acc_links_payload from accelerator_store collection
  struct get_acc_link_payload {
    template <typename acc_group_t, typename index_t>
    constexpr auto operator()(const acc_group_t& /*unused*/,
                              const index_t& index) const {
      using accel_t = typename acc_group_t::value_type;

      auto id{acc_links_payload::type_id::unknown};

      // Only convert grids
      if constexpr (detray::concepts::surface_grid<accel_t>) {
        id = io::detail::get_id<accel_t>();
      }

      return detail::basic_converter::to_payload(id, index);
    }
  };
};

}  // namespace detray::io
