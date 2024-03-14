// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts {

class CylinderVolumeStack : public Volume {
 public:
  enum class AttachmentStrategy { First, Second, Midpoint, Gap };

  CylinderVolumeStack(
      std::vector<std::shared_ptr<Volume>>& volumes, BinningValue direction,
      AttachmentStrategy strategy = AttachmentStrategy::Midpoint,
      const Logger& logger = Acts::getDummyLogger());

  // @TODO: Implement assignVolumeBounds

 private:
  static Volume createOuterVolume(std::vector<std::shared_ptr<Volume>>& volumes,
                                  BinningValue direction,
                                  AttachmentStrategy strategy,
                                  const Logger& logger);

  struct VolumeTuple;
  static std::shared_ptr<Volume> checkOverlapAndAttachInZ(
      VolumeTuple& a, VolumeTuple& b, const Transform3& groupTransform,
      AttachmentStrategy strategy, const Logger& logger);

  static void printVolumeSequence(const std::vector<VolumeTuple>& volumes,
                                  const Logger& logger,
                                  Acts::Logging::Level lvl);

  BinningValue m_direction;
  std::vector<std::shared_ptr<Volume>>& m_volumes;
};

std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::AttachmentStrategy strategy);

}  // namespace Acts
