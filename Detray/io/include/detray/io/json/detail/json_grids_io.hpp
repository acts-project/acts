// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/grid_axis.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/json/detail/json_algebra_io.hpp"
#include "detray/io/json/detail/json_common_io.hpp"
#include "detray/io/json/json.hpp"

// System include(s).
#include <array>
#include <optional>
#include <vector>

namespace detray::io {

inline void to_json(nlohmann::ordered_json& j, const grid_header_payload& h) {
  j["common"] = h.common;

  if (h.sub_header.has_value()) {
    const auto& grid_sub_header = h.sub_header.value();
    j["grid_count"] = grid_sub_header.n_grids;
  }
}

inline void from_json(const nlohmann::ordered_json& j, grid_header_payload& h) {
  h.common = j["common"];

  if (j.find("grid_count") != j.end()) {
    h.sub_header.emplace();
    auto& grid_sub_header = h.sub_header.value();
    grid_sub_header.n_grids = j["grid_count"];
  }
}

inline void to_json(nlohmann::ordered_json& j, const axis_payload& a) {
  j["label"] = static_cast<unsigned int>(a.label);
  j["bounds"] = static_cast<unsigned int>(a.bounds);
  j["binning"] = static_cast<unsigned int>(a.binning);
  j["bins"] = a.bins;
  j["edges"] = a.edges;
}

inline void from_json(const nlohmann::ordered_json& j, axis_payload& a) {
  a.binning = static_cast<axis::binning>(j["binning"]);
  a.bounds = static_cast<axis::bounds>(j["bounds"]);
  a.label = static_cast<axis::label>(j["label"]);
  a.bins = j["bins"];
  a.edges = j["edges"].get<std::vector<io::scalar>>();
}

template <typename content_t>
inline void to_json(nlohmann::ordered_json& j,
                    const grid_bin_payload<content_t>& g) {
  j["loc_index"] = g.loc_index;
  j["content"] = g.content;
}

template <typename content_t>
inline void from_json(const nlohmann::ordered_json& j,
                      grid_bin_payload<content_t>& g) {
  g.loc_index = j["loc_index"].get<std::vector<unsigned int>>();
  g.content = j["content"].get<std::vector<content_t>>();
}

template <typename content_t, typename grid_id_t>
inline void to_json(nlohmann::ordered_json& j,
                    const grid_payload<content_t, grid_id_t>& g) {
  j["owner_link"] = g.owner_link;
  j["grid_link"] = g.grid_link;

  nlohmann::ordered_json jaxes;
  for (const auto& a : g.axes) {
    jaxes.push_back(a);
  }
  j["axes"] = jaxes;

  nlohmann::ordered_json jbins;
  for (const auto& bin : g.bins) {
    jbins.push_back(bin);
  }
  j["bins"] = jbins;
}

template <typename content_t, typename grid_id_t>
inline void from_json(const nlohmann::ordered_json& j,
                      grid_payload<content_t, grid_id_t>& g) {
  g.owner_link = j["owner_link"];
  g.grid_link = j["grid_link"];

  nlohmann::ordered_json jaxes = j["axes"];
  for (auto jax : jaxes) {
    axis_payload a = jax;
    g.axes.push_back(std::move(a));
  }

  nlohmann::ordered_json jbins = j["bins"];
  for (auto jbin : jbins) {
    grid_bin_payload<content_t> b = jbin;
    g.bins.push_back(std::move(b));
  }
}

template <typename content_t, typename grid_id_t>
inline void to_json(nlohmann::ordered_json& j,
                    const detector_grids_payload<content_t, grid_id_t>& d) {
  if (!d.grids.empty()) {
    // Collection of volumes with their grid content
    nlohmann::ordered_json jgrids;

    for (const auto& [vol_idx, grs] : d.grids) {
      // Grids and volume they belong to
      nlohmann::ordered_json jgrid;
      // Array of grid payloads
      nlohmann::ordered_json grid_data;

      for (const auto& gr : grs) {
        grid_data.push_back(gr);
      }

      jgrid["volume_link"] = vol_idx;
      jgrid["grid_data"] = grid_data;
      jgrids.push_back(jgrid);
    }

    j["grids"] = jgrids;
  }
}

template <typename content_t, typename grid_id_t>
inline void from_json(const nlohmann::ordered_json& j,
                      detector_grids_payload<content_t, grid_id_t>& d) {
  if (j.find("grids") != j.end()) {
    for (auto& jgrids : j["grids"]) {
      const std::size_t vol_idx = jgrids["volume_link"];

      for (auto jgrid : jgrids["grid_data"]) {
        grid_payload<content_t, grid_id_t> grp = jgrid;
        d.grids[vol_idx].push_back(std::move(grp));
      }
    }
  }
}

}  // namespace detray::io
