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

// Project include(s)
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/json/detail/json_common_io.hpp"
#include "detray/io/json/json.hpp"

// System include(s)
#include <array>

/// @brief  The detray JSON I/O is written in such a way that it
/// can read/write ACTS files that are written with the Detray
/// JSON I/O extension
namespace detray::io {

inline void to_json(nlohmann::ordered_json& j,
                    const homogeneous_material_header_payload& h) {
  j["common"] = h.common;

  if (h.sub_header.has_value()) {
    const auto& mat_sub_header = h.sub_header.value();
    j["slab_count"] = mat_sub_header.n_slabs;
    j["rod_count"] = mat_sub_header.n_rods;
    j["slab_surface_count"] = mat_sub_header.n_slab_surfaces;
    j["rod_surface_count"] = mat_sub_header.n_rod_surfaces;
  }
}

inline void from_json(const nlohmann::ordered_json& j,
                      homogeneous_material_header_payload& h) {
  h.common = j["common"];

  if (j.find("slab_count") != j.end() && j.find("rod_count") != j.end()) {
    h.sub_header.emplace();
    auto& mat_sub_header = h.sub_header.value();
    mat_sub_header.n_slabs = j["slab_count"];
    mat_sub_header.n_rods = j["rod_count"];
    // Default surface counts to 0 to keep older files compatible
    mat_sub_header.n_slab_surfaces =
        j.value<std::size_t>("slab_surface_count", 0);
    mat_sub_header.n_rod_surfaces =
        j.value<std::size_t>("rod_surface_count", 0);
  }
}

inline void to_json(nlohmann::ordered_json& j,
                    const material_param_payload& m) {
  j["params"] = m.params;
}

inline void from_json(const nlohmann::ordered_json& j,
                      material_param_payload& m) {
  m.params = j["params"].get<std::array<io::scalar, 7>>();
}

inline void to_json(nlohmann::ordered_json& j,
                    const surface_material_payload& m) {
  j["type"] = m.type;
  j["surface_idx"] = m.surface;
  j["thickness"] = m.thickness;
  j["material"] = m.mat;
  if (m.index_in_coll.has_value()) {
    j["index_in_coll"] = m.index_in_coll.value();
  }
}

inline void from_json(const nlohmann::ordered_json& j,
                      surface_material_payload& m) {
  m.type = j["type"];
  m.surface = j["surface_idx"];
  m.thickness = j["thickness"];
  m.mat = j["material"];
  if (j.find("index_in_coll") != j.end()) {
    m.index_in_coll = j["index_in_coll"];
  }
}

inline void to_json(nlohmann::ordered_json& j,
                    const material_volume_payload& mv) {
  j["volume_link"] = mv.volume_link;

  if (!mv.surface_mat.empty()) {
    nlohmann::ordered_json jmats;
    for (const auto& m : mv.surface_mat) {
      jmats.push_back(m);
    }
    j["surface_material"] = jmats;
  }
}

inline void from_json(const nlohmann::ordered_json& j,
                      material_volume_payload& mv) {
  mv.volume_link = j["volume_link"];

  if (j.find("surface_material") != j.end()) {
    for (auto jmats : j["surface_material"]) {
      surface_material_payload mslp = jmats;
      mv.surface_mat.push_back(mslp);
    }
  } else if (j.find("material_slabs") != j.end()) {
    // TODO: Remove once all files are fixed
    for (auto jmats : j["material_slabs"]) {
      surface_material_payload mslp = jmats;
      mv.surface_mat.push_back(mslp);
    }
  }
}

inline void to_json(nlohmann::ordered_json& j,
                    const detector_homogeneous_material_payload& d) {
  if (!d.volumes.empty()) {
    nlohmann::ordered_json jmats;
    for (const auto& m : d.volumes) {
      jmats.push_back(m);
    }
    j["volumes"] = jmats;
  }
}

inline void from_json(const nlohmann::ordered_json& j,
                      detector_homogeneous_material_payload& d) {
  if (j.find("volumes") != j.end()) {
    for (auto jvolume : j["volumes"]) {
      d.volumes.push_back(jvolume);
    }
  }
}

}  // namespace detray::io
