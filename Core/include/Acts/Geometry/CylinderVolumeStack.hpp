// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
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
  Volume createOuterVolume(std::vector<std::shared_ptr<Volume>>& volumes,
                           BinningValue direction, AttachmentStrategy strategy,
                           const Logger& logger);

  struct VolumeTuple;

  static void printVolumeSequence(const std::vector<VolumeTuple>& volumes,
                                  const Logger& logger,
                                  Acts::Logging::Level lvl);

  static void overlapPrint(const VolumeTuple& a, const VolumeTuple& b,
                           const Logger& logger);

  static void checkVolumeAlignment(const std::vector<VolumeTuple>& volumes,
                                   const Logger& logger);

  std::vector<VolumeTuple> checkOverlapAndAttachInZ(
      std::vector<VolumeTuple>& volumes, AttachmentStrategy strategy,
      const Logger& logger);

  std::pair<ActsScalar, ActsScalar> synchronizeRBounds(
      std::vector<VolumeTuple>& volumes, const Logger& logger);

  std::vector<VolumeTuple> checkOverlapAndAttachInR(
      std::vector<VolumeTuple>& volumes, AttachmentStrategy strategy,
      const Logger& logger);

  std::pair<ActsScalar, ActsScalar> synchronizeZBounds(
      std::vector<VolumeTuple>& volumes, const Logger& logger);

  BinningValue m_direction{};
  Transform3 m_groupTransform{};
  std::vector<std::shared_ptr<Volume>>& m_volumes;
};

std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::AttachmentStrategy strategy);

}  // namespace Acts
