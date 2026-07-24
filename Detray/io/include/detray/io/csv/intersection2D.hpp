// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/io/csv/dfe.hpp"
#include "detray/io/utils/create_path.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/utils/ranges.hpp"

// System include(s).
#include <cstdint>
#include <filesystem>

namespace detray::io::csv {

/// Type to read the data of a line-surface intersection
struct intersection2D {
  unsigned int track_id = 0;
  std::uint64_t identifier = 0ul;
  unsigned int type = 0u;
  unsigned int transform_index = 0u;
  unsigned int mask_id = 0u;
  unsigned int mask_index = 0u;
  unsigned int n_masks = 1u;
  unsigned int material_id = 0u;
  unsigned int material_index = 0u;
  double loc_0 = 0.;
  double loc_1 = 0.;
  double path = 0.;
  unsigned int volume_link = 0u;
  int direction = 0;
  int status = 0;

  DETRAY_DFE_NAMEDTUPLE(intersection2D, track_id, identifier, type,
                        transform_index, mask_id, mask_index, material_id,
                        material_index, loc_0, loc_1, path, volume_link,
                        direction, status);
};

/// Read intersections from csv file
/// @returns vector of intersections
template <typename detector_t>
inline auto read_intersection2D(const std::string &file_name) {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using surface_t = typename detector_t::surface_type;
  using nav_link_t = typename surface_t::navigation_link;
  using mask_link_t = typename surface_t::mask_link;
  using material_link_t = typename surface_t::material_link;
  using mask_id_t = typename detector_t::masks::id;
  using material_id_t = typename detector_t::material::id;

  using intersection_t =
      detray::intersection2D<surface_t, algebra_t, intersection::contains_pos>;

  dfe::NamedTupleCsvReader<io::csv::intersection2D> inters_reader(file_name);

  io::csv::intersection2D inters_data{};
  std::vector<std::vector<intersection_t>> intersections_per_track;

  // Read the data
  while (inters_reader.read(inters_data)) {
    // Add new intersection to correct track
    dindex trk_index{inters_data.track_id};
    while (intersections_per_track.size() <= trk_index) {
      intersections_per_track.push_back({});
    }

    // Read the intersection
    intersection_t inters{};

    mask_link_t mask_link{};
    mask_link.set_id(static_cast<mask_id_t>(inters_data.mask_id));
    if constexpr (detray::concepts::interval<
                      typename mask_link_t::index_type>) {
      mask_link.set_index({inters_data.mask_index, inters_data.n_masks});
    } else {
      mask_link.set_index(inters_data.mask_index);
    }

    material_link_t material_link{
        static_cast<material_id_t>(inters_data.material_id),
        inters_data.material_index};
    inters.set_surface({inters_data.transform_index, mask_link, material_link,
                        dindex_invalid, surface_id::e_unknown});
    inters.surface().set_identifier(
        geometry::identifier{inters_data.identifier});
    inters.set_local({static_cast<scalar_t>(inters_data.loc_0),
                      static_cast<scalar_t>(inters_data.loc_1), 0.f});
    inters.set_path(static_cast<scalar_t>(inters_data.path));
    inters.set_volume_link(static_cast<nav_link_t>(inters_data.volume_link));
    inters.set_direction(static_cast<bool>(inters_data.direction));
    inters.set_status(static_cast<intersection::status>(inters_data.status));

    // Add to collection
    intersections_per_track[trk_index].push_back(inters);
  }

  // Check the result
  if (intersections_per_track.empty()) {
    throw std::invalid_argument(
        "ERROR: csv reader: Failed to read intersection data");
  }

  return intersections_per_track;
}

/// Write intersections to csv file
template <typename intersection_t>
inline void write_intersection2D(
    const std::string &file_name,
    const std::vector<std::vector<intersection_t>> &intersections_per_track,
    const bool replace = true) {
  using sf_desc_t = decltype(intersections_per_track.at(0).at(0).surface());
  using mask_link_t = typename sf_desc_t::mask_link;

  // Don't write over existing data
  std::string inters_file_name{file_name};
  if (!replace && io::file_exists(file_name)) {
    inters_file_name = io::alt_file_name(file_name);
  } else {
    // Make sure the output directories exit
    io::create_path(std::filesystem::path{inters_file_name}.parent_path());
  }

  dfe::NamedTupleCsvWriter<io::csv::intersection2D> inters_writer(
      inters_file_name);

  for (const auto &[track_idx, intersections] :
       detray::views::enumerate(intersections_per_track)) {
    // Skip empty traces
    if (intersections.empty()) {
      continue;
    }

    for (const auto &inters : intersections) {
      io::csv::intersection2D inters_data{};

      inters_data.track_id = track_idx;
      inters_data.identifier = inters.surface().identifier().value();
      inters_data.type =
          static_cast<unsigned int>(inters.surface().identifier().id());
      inters_data.transform_index = inters.surface().transform();
      inters_data.mask_id =
          static_cast<unsigned int>(inters.surface().mask().id());
      const auto mask_index = inters.surface().mask().index();
      if constexpr (detray::concepts::interval<
                        typename mask_link_t::index_type>) {
        inters_data.mask_index = static_cast<unsigned int>(mask_index.lower());
        inters_data.n_masks = static_cast<unsigned int>(mask_index.size());
      } else {
        inters_data.mask_index = static_cast<unsigned int>(mask_index);
        inters_data.n_masks = 1u;
      }
      inters_data.material_id =
          static_cast<unsigned int>(inters.surface().material().id());
      inters_data.material_index = inters.surface().material().index();
      inters_data.loc_0 = inters.local()[0];
      inters_data.loc_1 = inters.local()[1];
      inters_data.path = inters.path();
      inters_data.volume_link = inters.volume_link();
      inters_data.direction = static_cast<int>(inters.is_along());
      inters_data.status = static_cast<int>(inters.status());

      inters_writer.append(inters_data);
    }
  }
}

}  // namespace detray::io::csv
