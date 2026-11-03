// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/utils/io_metadata.hpp"

// System include(s)
#include <string_view>

// Convert basic information like links and header data to and from payloads
namespace detray::io::detail::basic_converter {

/// @returns a link from its io payload @param link_data
inline dindex from_payload(const single_link_payload& link_data) {
  return static_cast<dindex>(link_data.link);
}

/// Convert a link @param idx into its io payload
inline single_link_payload to_payload(const std::size_t idx) {
  single_link_payload link_data;
  link_data.link = idx;

  return link_data;
}

/// Convert a typed link with a type id @param id and index @param idx into its
/// io payload
template <typename type_id>
inline typed_link_payload<type_id> to_payload(const type_id id,
                                              const std::size_t idx) {
  typed_link_payload<type_id> link_data;

  link_data.type = id;
  link_data.index = idx;

  return link_data;
}

/// Convert the common header information using the detector name
/// @param det_name and the file tag @param tag that describes the data file
/// content
template <detray::concepts::detector detector_t>
inline tagged_header_payload to_header_payload(const std::string_view det_name,
                                               const std::string_view tag,
                                               const std::string_view source) {
  tagged_header_payload header_data;

  header_data.version = io::detail::get_detray_version();
  header_data.date = io::detail::get_current_date();
  header_data.metadata =
      detray::types::get_name<typename detector_t::metadata>();
  header_data.source = source;
  header_data.detector = det_name;
  header_data.tag = tag;

  return header_data;
}

/// @returns the geometry header data, generated from the geometry payload
/// @param geo_data
inline geometry_header_payload to_header_payload(
    const detector_geometry_payload& geo_data) {
  geometry_header_payload header_data;

  header_data.n_volumes = geo_data.volumes.size();
  for (const auto& vol_data : geo_data.volumes) {
    header_data.n_surfaces += vol_data.surfaces.size();
  }

  return header_data;
}

/// @returns the homogeneous material header data, generated from the
/// homogeneous material payload @param mat_data
inline homogeneous_material_header_payload to_header_payload(
    const detector_homogeneous_material_payload& mat_data) {
  using material_id_t = surface_material_payload::mat_type;

  homogeneous_material_header_payload header_data;

  // Every surface material payload belongs to exactly one surface
  for (const auto& mv_data : mat_data.volumes) {
    for (const auto& mslp : mv_data.surface_mat) {
      if (mslp.type == material_id_t::slab) {
        ++header_data.n_slabs;
      } else if (mslp.type == material_id_t::rod) {
        ++header_data.n_rods;
      }
    }
  }
  header_data.n_slab_surfaces = header_data.n_slabs;
  header_data.n_rod_surfaces = header_data.n_rods;

  return header_data;
}

/// @returns the grid header data, generated from the grid collection payload
/// @param grids_data
template <typename content_t, typename grid_id_t>
inline grid_header_payload to_header_payload(
    const detector_grids_payload<content_t, grid_id_t>& grids_data) {
  grid_header_payload header_data;

  for (const auto& [vol_idx, grids] : grids_data.grids) {
    header_data.n_grids += grids.size();
  }

  return header_data;
}

/// Fill the header of the detector payload @param payload: Generate the
/// common information and the sub-headers from the contained payload data
///
/// @param source the original data source of the detector. Does not
///               overwrite a known source, if it is 'unknown'
inline void fill_header_payload(detector_payload& payload,
                                const std::string_view source) {
  header_payload& header = payload.header;

  header.version = io::detail::get_detray_version();
  header.date = io::detail::get_current_date();
  if (header.metadata.empty()) {
    header.metadata = "unknown";
  }
  if (header.source.empty() || source != "unknown") {
    header.source = source;
  }

  header.geometry = to_header_payload(payload.geometry);
  if (payload.volume_grids.has_value()) {
    header.volume_grids = to_header_payload(*payload.volume_grids);
  }
  if (payload.surface_grids.has_value()) {
    header.surface_grids = to_header_payload(*payload.surface_grids);
  }
  if (payload.homogeneous_material.has_value()) {
    header.homogeneous_material =
        to_header_payload(*payload.homogeneous_material);
  }
  if (payload.material_maps.has_value()) {
    header.material_maps = to_header_payload(*payload.material_maps);
  }
}

}  // namespace detray::io::detail::basic_converter
