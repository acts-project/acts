// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/io/csv/dfe.hpp"
#include "detray/io/utils/create_path.hpp"
#include "detray/tracks/free_track_parameters.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <filesystem>

namespace detray::io::csv {

/// Type to read the data of free track parameters
struct free_track_parameters {
  unsigned int track_id = 0;
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;
  double px = 0.;
  double py = 0.;
  double pz = 0.;
  double q = 0.;

  DETRAY_DFE_NAMEDTUPLE(free_track_parameters, track_id, x, y, z, t, px, py, pz,
                        q);
};

/// Type to read the data of bound track parameters
/// @TODO: Implement readers and writers
struct bound_track_parameters {
  unsigned int track_id = 0;
  double l0 = 0.;
  double l1 = 0.;
  double phi = 0.;
  double theta = 0.;
  double qop = 0.;
  double t = 0.;

  DETRAY_DFE_NAMEDTUPLE(bound_track_parameters, track_id, l0, l1, phi, theta,
                        qop, t);
};

/// Read free track parameters from csv file
/// @returns vector of free track parameters
template <typename detector_t>
inline auto read_free_track_params(const std::string &file_name) {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;
  using track_t = detray::free_track_parameters<algebra_t>;

  dfe::NamedTupleCsvReader<io::csv::free_track_parameters> track_param_reader(
      file_name);

  io::csv::free_track_parameters track_param_data{};
  std::vector<std::vector<std::pair<scalar_t, track_t>>> track_params_per_track;

  while (track_param_reader.read(track_param_data)) {
    // Add new track state to correct track
    dindex trk_index{track_param_data.track_id};
    while (track_params_per_track.size() <= trk_index) {
      track_params_per_track.emplace_back();
    }

    point3_t pos{static_cast<scalar_t>(track_param_data.x),
                 static_cast<scalar_t>(track_param_data.y),
                 static_cast<scalar_t>(track_param_data.z)};
    vector3_t p{static_cast<scalar_t>(track_param_data.px),
                static_cast<scalar_t>(track_param_data.py),
                static_cast<scalar_t>(track_param_data.pz)};

    track_t track_param{pos, static_cast<scalar_t>(track_param_data.t), p,
                        static_cast<scalar_t>(track_param_data.q)};

    // Add to collection
    track_params_per_track[trk_index].emplace_back(track_param_data.q,
                                                   track_param);
  }

  // Check the result
  if (track_params_per_track.empty()) {
    throw std::invalid_argument(
        "ERROR: csv reader: Failed to read free track parameters data");
  }

  return track_params_per_track;
}

/// Write free track parameters to csv file
template <detray::concepts::scalar scalar_t, typename track_t>
inline void write_free_track_params(
    const std::string &file_name,
    const std::vector<std::vector<std::pair<scalar_t, track_t>>>
        &track_params_per_track,
    const bool replace = true) {
  // Don't write over existing data
  std::string trk_file_name{file_name};
  if (!replace && io::file_exists(file_name)) {
    trk_file_name = io::alt_file_name(file_name);
  } else {
    // Make sure the output directories exit
    io::create_path(std::filesystem::path{trk_file_name}.parent_path());
  }

  dfe::NamedTupleCsvWriter<io::csv::free_track_parameters> track_param_writer(
      trk_file_name);

  for (const auto &[track_idx, track_params] :
       detray::views::enumerate(track_params_per_track)) {
    // Skip empty traces
    if (track_params.empty()) {
      continue;
    }

    for (const auto &[charge, track_param] : track_params) {
      const auto &glob_pos = track_param.pos();
      // Momentum may not be retrievable for straight-line tracks
      const auto &p{charge != 0.f ? track_param.mom(charge)
                                  : track_param.dir()};

      io::csv::free_track_parameters track_param_data{};
      track_param_data.track_id = track_idx;
      track_param_data.x = glob_pos[0];
      track_param_data.y = glob_pos[1];
      track_param_data.z = glob_pos[2];
      track_param_data.t = track_param.time();
      track_param_data.px = p[0];
      track_param_data.py = p[1];
      track_param_data.pz = p[2];
      track_param_data.q = charge;

      track_param_writer.append(track_param_data);
    }
  }
}

}  // namespace detray::io::csv
