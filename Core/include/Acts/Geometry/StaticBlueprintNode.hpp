// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

class StaticBlueprintNode : public BlueprintNode {
 public:
  StaticBlueprintNode(const std::string& name,
                      std::unique_ptr<TrackingVolume> volume);

  Volume& build() override;

  void connect(TrackingVolume& parent) override;

  // This connects averything to the static volume contained here
  void connect();

  std::unique_ptr<TrackingVolume> releaseVolume() {
    return std::move(m_volume);
  }

 private:
  std::unique_ptr<TrackingVolume> m_volume;
};

}  // namespace Acts
