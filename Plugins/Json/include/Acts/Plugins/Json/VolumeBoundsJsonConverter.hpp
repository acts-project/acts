// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {
class VolumeBounds;

const static std::vector<std::string> volumeBoundTypes = {
    "Cone",     "Cuboid",        "CutoutCylinder",
    "Cylinder", "GenericCuboid", "Trapezoid"};

void to_json(nlohmann::json& j, const VolumeBounds& bounds);

void to_json(nlohmann::json& j, const GenericCuboidVolumeBounds& bounds);

/// Conversion to surfaceBounds from json
///
/// The type is given as a template argument in order to be able
/// to construct the correct fitting types for surfaces.
///
/// @param j the read-in json object
///
/// @return a shared_ptr to a surface object for type polymorphism
template <typename bounds_t>
std::unique_ptr<bounds_t> volumeBoundsFromJson(const nlohmann::json& j) {
  const size_t kValues = bounds_t::BoundValues::eSize;
  std::array<ActsScalar, kValues> bValues{};
  std::vector<ActsScalar> bVector = j["values"];
  std::copy_n(bVector.begin(), kValues, bValues.begin());
  return std::make_unique<bounds_t>(bValues);
}

inline std::unique_ptr<GenericCuboidVolumeBounds> genericVolumeBoundsFromJson(
    const nlohmann::json& j) {
  auto json_vertices = j["values"];
  std::array<Vector3, 8> vertices;
  for (size_t i = 0; i < 8; i++) {
    vertices[i] << json_vertices[i][0], json_vertices[i][1],
        json_vertices[i][2];
  }
  return std::make_unique<GenericCuboidVolumeBounds>(vertices);
}

/// Conversion to surfaceBounds from json
/// @param j the read-in json object
///
/// @return a shared_ptr to a surface object for type polymorphism
std::unique_ptr<VolumeBounds> unqiueVolumeBoundsFromJson(
    const nlohmann::json& j);

}  // namespace Acts
