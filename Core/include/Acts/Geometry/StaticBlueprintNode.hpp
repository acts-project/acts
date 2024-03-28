// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

class StaticBlueprintNode : public BlueprintNode {
 public:
  StaticBlueprintNode(std::unique_ptr<TrackingVolume> volume);

  const std::string& name() const override;

  Volume& build(const Logger& logger = Acts::getDummyLogger()) override;

  void connect(TrackingVolume& parent,
               const Logger& logger = Acts::getDummyLogger()) override;

  void visualize(IVisualization3D& vis,
                 const GeometryContext& gctx) const override;

  // This connects averything to the static volume contained here
  void connect(const Logger& logger = Acts::getDummyLogger());

  std::unique_ptr<TrackingVolume> releaseVolume() {
    return std::move(m_volume);
  }

  // protected:
  // void addToGraphviz(std::ostream& os) const override;

 private:
  std::unique_ptr<TrackingVolume> m_volume;
};

}  // namespace Acts
