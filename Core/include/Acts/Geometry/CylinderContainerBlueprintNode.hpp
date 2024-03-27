// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"
#include "Acts/Utilities/BinningType.hpp"
namespace Acts {

class CylinderContainerBlueprintNode : public BlueprintNode {
 public:
  CylinderContainerBlueprintNode(
      const std::string& name, BinningValue direction,
      CylinderVolumeStack::AttachmentStrategy attachmentStrategy =
          CylinderVolumeStack::AttachmentStrategy::Midpoint,
      CylinderVolumeStack::ResizeStrategy resizeStrategy =
          CylinderVolumeStack::ResizeStrategy::Expand);

  Volume& build() override;

  void connect(TrackingVolume& parent) override;

 private:
  BinningValue m_direction;
  CylinderVolumeStack::AttachmentStrategy m_attachmentStrategy;
  CylinderVolumeStack::ResizeStrategy m_resizeStrategy;
  // CylinderVolumeStack m_stack;
};

}  // namespace Acts
