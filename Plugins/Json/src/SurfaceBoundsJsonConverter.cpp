// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/SurfaceBoundsJsonConverter.hpp"

#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "ActsPlugins/Json/DetrayJsonHelper.hpp"

void Acts::to_json(nlohmann::json& j, const Acts::SurfaceBounds& bounds) {
  j["type"] = bounds.type();
  j["values"] = bounds.values();
}

nlohmann::json Acts::SurfaceBoundsJsonConverter::toJson(
    const Acts::SurfaceBounds& bounds) {
  return nlohmann::json(bounds);
}

nlohmann::json Acts::SurfaceBoundsJsonConverter::toJsonDetray(
    const Acts::SurfaceBounds& bounds, bool portal) {
  nlohmann::json jMask;
  auto [shape, boundaries] = DetrayJsonHelper::maskFromBounds(bounds, portal);
  jMask["shape"] = shape;
  jMask["boundaries"] = boundaries;
  return jMask;
}
