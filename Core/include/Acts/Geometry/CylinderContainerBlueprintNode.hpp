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
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>

namespace Acts {

class CylinderContainerBlueprintNode final : public BlueprintNode {
 public:
  CylinderContainerBlueprintNode(
      const std::string& name, BinningValue direction,
      CylinderVolumeStack::AttachmentStrategy attachmentStrategy =
          CylinderVolumeStack::AttachmentStrategy::Midpoint,
      CylinderVolumeStack::ResizeStrategy resizeStrategy =
          CylinderVolumeStack::ResizeStrategy::Expand);

  const std::string& name() const override;

  void setName(const std::string& name) { m_name = name; }

  Volume& build(const Logger& logger = Acts::getDummyLogger()) override;

  CylinderStackPortalShell& connect(
      const GeometryContext& gctx,
      const Logger& logger = Acts::getDummyLogger()) override;

  void finalize(TrackingVolume& parent, const Logger& logger) override;

  CylinderContainerBlueprintNode& setDirection(BinningValue direction);
  CylinderContainerBlueprintNode& setAttachmentStrategy(
      CylinderVolumeStack::AttachmentStrategy attachmentStrategy);
  CylinderContainerBlueprintNode& setResizeStrategy(
      CylinderVolumeStack::ResizeStrategy resizeStrategy);

  BinningValue direction() const { return m_direction; }
  CylinderVolumeStack::AttachmentStrategy attachmentStrategy() const {
    return m_attachmentStrategy;
  }
  CylinderVolumeStack::ResizeStrategy resizeStrategy() const {
    return m_resizeStrategy;
  }

 private:
  void addToGraphviz(std::ostream& os) const override;

  bool isGapVolume(const Volume& volume) const;

  std::string m_name;

  BinningValue m_direction;
  CylinderVolumeStack::AttachmentStrategy m_attachmentStrategy;
  CylinderVolumeStack::ResizeStrategy m_resizeStrategy;

  // Is only initialized during `build`
  std::vector<Volume*> m_childVolumes;
  std::unique_ptr<CylinderVolumeStack> m_stack{nullptr};
  std::map<const Volume*, BlueprintNode*> m_volumeToNode;

  std::vector<std::pair<std::unique_ptr<SingleCylinderPortalShell>,
                        std::unique_ptr<TrackingVolume>>>
      m_gaps;

  std::unique_ptr<CylinderStackPortalShell> m_shell{nullptr};
};

}  // namespace Acts
