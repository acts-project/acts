// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

namespace Acts {

class CylinderVolumeStack : public Volume {
 public:
  enum class AttachmentStrategy { First, Second, Midpoint, Gap };
  enum class ResizeStrategy { Expand, Gap };

  CylinderVolumeStack(
      std::vector<Volume*>& volumes, BinningValue direction,
      AttachmentStrategy strategy = AttachmentStrategy::Midpoint,
      ResizeStrategy resizeStrategy = ResizeStrategy::Expand,
      const Logger& logger = Acts::getDummyLogger());

  void assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds) override;

  void assignVolumeBounds(std::shared_ptr<VolumeBounds> volbounds,
                          const Logger& logger);

 private:
  static Volume initialVolume(const std::vector<Volume*>& volumes);

  void initializeOuterVolume(BinningValue direction,
                             AttachmentStrategy strategy, const Logger& logger);

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

  static void checkNoPhiOrBevel(const CylinderVolumeBounds& bounds,
                                const Logger& logger);

  std::shared_ptr<Volume> addGapVolume(Transform3 transform,
                                       std::shared_ptr<VolumeBounds> bounds);

  std::vector<std::shared_ptr<Volume>>& gaps();

  BinningValue m_direction{};
  ResizeStrategy m_resizeStrategy{};
  Transform3 m_groupTransform{};
  std::vector<std::shared_ptr<Volume>> m_gaps{};
  std::vector<Volume*>& m_volumes;
};

std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::AttachmentStrategy strategy);
std::ostream& operator<<(std::ostream& os,
                         CylinderVolumeStack::ResizeStrategy strategy);

}  // namespace Acts
