// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "materialPlotHelper.hpp"

#include <iomanip>
#include <ostream>
#include <string>
#include <unordered_map>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

/// Information on a given surface.

struct sinfo {
  std::string name;
  std::string idname;
  std::string id;
  int type;
  float pos;
  float range_min;
  float range_max;
};

std::ostream& Acts::operator<<(std::ostream& os, Acts::GeometryIdentifier id) {
  os << "[ " << std::setw(3) << id.volume();
  os << " | " << std::setw(3) << id.boundary();
  os << " | " << std::setw(3) << id.layer();
  os << " | " << std::setw(3) << id.approach();
  os << " | " << std::setw(4) << id.sensitive() << " ]";
  return os;
}

std::unordered_map<std::uint64_t, json> load_geometry_file(
    std::string geometry_file) {
  json geom;
  {
    std::ifstream gj(geometry_file);
    if (!gj.good()) {
      std::cerr << "WARNING: " << geometry_file << " not found." << std::endl;
    } else {
      try {
        gj >> geom;
      } catch (...) {
        std::cerr << "WARNING: Failed to parse " << geometry_file << "."
                  << std::endl;
      }
    }
  }
  std::unordered_map<std::uint64_t, json> surface_bounds;
  const auto& entries = geom["Surfaces"]["entries"];
  for (const auto& entry : entries) {
    // Handle both old (uint64_t) and new (object with component fields) formats
    std::uint64_t gid;
    if (entry["value"]["geo_id"].is_number()) {
      // Specification: geo_id is a simple uint64_t (
      gid = entry["value"]["geo_id"].get<std::uint64_t>();
    } else if (entry["value"]["geo_id"].is_object()) {
      // Specification: geo_id is an object with component fields
      // Reconstruct the uint64_t from individual components
      // Layout: volume (bits 56-63), boundary (bits 48-55), layer (bits 36-47),
      //         approach (bits 28-35), sensitive (bits 0-27)
      std::uint64_t gid_value = 0;
      const auto& geo_id_obj = entry["value"]["geo_id"];

      if (geo_id_obj.contains("volume") && !geo_id_obj["volume"].is_null()) {
        uint32_t volume = static_cast<uint32_t>(geo_id_obj["volume"].get<double>());
        gid_value |= (static_cast<std::uint64_t>(volume) & 0xFF) << 56;
      }
      if (geo_id_obj.contains("boundary") && !geo_id_obj["boundary"].is_null()) {
        uint32_t boundary = static_cast<uint32_t>(geo_id_obj["boundary"].get<double>());
        gid_value |= (static_cast<std::uint64_t>(boundary) & 0xFF) << 48;
      }
      if (geo_id_obj.contains("layer") && !geo_id_obj["layer"].is_null()) {
        uint32_t layer = static_cast<uint32_t>(geo_id_obj["layer"].get<double>());
        gid_value |= (static_cast<std::uint64_t>(layer) & 0xFFF) << 36;
      }
      if (geo_id_obj.contains("approach") && !geo_id_obj["approach"].is_null()) {
        uint32_t approach = static_cast<uint32_t>(geo_id_obj["approach"].get<double>());
        gid_value |= (static_cast<std::uint64_t>(approach) & 0xFF) << 28;
      }
      if (geo_id_obj.contains("sensitive") && !geo_id_obj["sensitive"].is_null()) {
        uint32_t sensitive = static_cast<uint32_t>(geo_id_obj["sensitive"].get<double>());
        gid_value |= (static_cast<std::uint64_t>(sensitive) & 0x0FFFFFFF);
      }
      gid = gid_value;
    } else {
      std::cerr << "WARNING: Unknown geo_id format in geometry file." << std::endl;
      continue;
    }
    surface_bounds[gid] = entry["value"]["bounds"];
  }
  return surface_bounds;
}

/// Initialise the information on each surface.

void Initialise_info(sinfo& surface_info,
                     const std::map<std::string, std::string>& surfaceName,
                     const std::uint64_t& id, const int& type, const float& pos,
                     const float& range_min, const float& range_max) {
  Acts::GeometryIdentifier ID(id);
  std::ostringstream layerID;
  layerID << ID;
  std::string surface_id = layerID.str();

  std::string Id_temp = surface_id;
  std::string delimiter = " | ";
  std::size_t del_pos = 0;
  std::vector<std::string> Ids;
  while ((del_pos = Id_temp.find(delimiter)) != std::string::npos) {
    Ids.push_back(Id_temp.substr(0, del_pos));
    Id_temp.erase(0, del_pos + delimiter.length());
  }
  Ids.push_back(Id_temp);

  for (int tag = 0; tag < 5; tag++) {
    Ids[tag].erase(std::remove(Ids[tag].begin(), Ids[tag].end(), ' '),
                   Ids[tag].end());
    Ids[tag].erase(std::remove(Ids[tag].begin(), Ids[tag].end(), '['),
                   Ids[tag].end());
    Ids[tag].erase(std::remove(Ids[tag].begin(), Ids[tag].end(), ']'),
                   Ids[tag].end());
  }

  surface_info.idname = "v" + Ids[0] + "_b" + Ids[1] + "_l" + Ids[2] + "_a" +
                        Ids[3] + "_s" + Ids[4];
  surface_info.type = type;

  if (surfaceName.contains(surface_id)) {
    surface_info.name = surfaceName.at(surface_id);
  } else {
    surface_info.name = "";
  }

  surface_info.id = surface_id;
  surface_info.pos = pos;
  surface_info.range_min = range_min;
  surface_info.range_max = range_max;
}
