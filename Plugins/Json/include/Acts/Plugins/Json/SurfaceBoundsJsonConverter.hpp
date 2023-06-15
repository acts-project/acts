// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {
class SurfaceBounds;

static std::vector<std::string> boundTypes = {
    "ConeBounds",          "CylinderBounds",      "DiamondBounds",
    "RadialBounds",        "EllipseBounds",       "LineBounds",
    "RectangleBounds",     "TrapezoidBounds",     "TriangleBounds",
    "DiscTrapezoidBounds", "ConvexPolygonBounds", "AnnulusBounds",
    "OtherBounds"};

void to_json(nlohmann::json& j, const SurfaceBounds& bounds);

/// Converstion to surfaceBounds from json
///
/// The type is given as a template argument in order to be able
/// to construct the correct fitting types for surfaces.
///
/// @param j the read-in json object
///
/// @return a shared_ptr to a surface object for type polymorphism
template <typename bounds_t>
std::shared_ptr<const bounds_t> surfaceBoundsFromJson(const nlohmann::json& j) {
  const size_t kValues = bounds_t::BoundValues::eSize;
  std::array<ActsScalar, kValues> bValues{};
  std::vector<ActsScalar> bVector = j["values"];
  std::copy_n(bVector.begin(), kValues, bValues.begin());
  return std::make_shared<const bounds_t>(bValues);
}

}  // namespace Acts