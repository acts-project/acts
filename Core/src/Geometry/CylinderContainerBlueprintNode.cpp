// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderContainerBlueprintNode.hpp"

namespace Acts {

CylinderContainerBlueprintNode::CylinderContainerBlueprintNode(
    const std::string& name, BinningValue direction,
    CylinderVolumeStack::AttachmentStrategy attachmentStrategy,
    CylinderVolumeStack::ResizeStrategy resizeStrategy)
    : BlueprintNode(name),
      m_direction(direction),
      m_attachmentStrategy(attachmentStrategy),
      m_resizeStrategy(resizeStrategy) {}

Volume& CylinderContainerBlueprintNode::build() {}

void CylinderContainerBlueprintNode::connect(TrackingVolume& parent) {}

}  // namespace Acts
