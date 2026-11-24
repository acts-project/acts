// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/VolumeJsonConverter.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "ActsPlugins/Json/GeometryJsonKeys.hpp"
#include "ActsPlugins/Json/MaterialJsonConverter.hpp"

void Acts::to_json(
    nlohmann::json& j,
    const std::pair<const Acts::TrackingVolume*,
                    std::shared_ptr<const Acts::IVolumeMaterial>>& volume) {
  j[Acts::jsonKey().namekey] = volume.first->volumeName();
  // Volume material has no replacement in the surface-only v1 format.
  ACTS_PUSH_IGNORE_DEPRECATED()
  to_json(j, volume.second.get());
  ACTS_POP_IGNORE_DEPRECATED()
}

void Acts::to_json(nlohmann::json& j, const Acts::TrackingVolume& volume) {
  j[Acts::jsonKey().namekey] = volume.volumeName();
  if (volume.hasMaterial()) {
    // Volume material has no replacement in the surface-only v1 format.
    ACTS_PUSH_IGNORE_DEPRECATED()
    to_json(j, volume.volumeMaterial());
    ACTS_POP_IGNORE_DEPRECATED()
  }
  return;
}
