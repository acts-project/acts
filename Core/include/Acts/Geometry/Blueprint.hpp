// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

namespace Acts {

class TrackingVolume;

class Blueprint {
  std::unique_ptr<TrackingVolume> build();

  // void setStaticVolume(std::unique_ptr<TrackingVolume> volume);

  // void setCylinderContainer(
  // const std::function<std::unique_ptr<BlueprintNode>(
  // std::unique_ptr<CylinderContainerBlueprintNode> cylinder)>& factory);
};

}  // namespace Acts
