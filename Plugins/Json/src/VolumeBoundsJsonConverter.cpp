// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"

#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::VolumeBounds& bounds) {
  j["type"] = volumeBoundTypes[bounds.type()];
  j["values"] = bounds.values();
}

std::unique_ptr<Acts::VolumeBounds> Acts::unqiueVolumeBoundsFromJson(
    const nlohmann::json& j) {
  const std::string type = j["type"];

  if (type == "Cone") {
    return volumeBoundsFromJson<ConeVolumeBounds>(j);
  } else if (type == "Cuboid") {
    return volumeBoundsFromJson<CuboidVolumeBounds>(j);
  } else if (type == "CutoutCylinder") {
    return volumeBoundsFromJson<CutoutCylinderVolumeBounds>(j);
  } else if (type == "Cylinder") {
    return volumeBoundsFromJson<CylinderVolumeBounds>(j);
  } else if (type == "Trapezoid") {
    return volumeBoundsFromJson<TrapezoidVolumeBounds>(j);
  }
  // case "GenericCuboid": {
  //   return volumeBoundsFromJson<GenericCuboidVolumeBounds>(j);
  // };
  return nullptr;
}
