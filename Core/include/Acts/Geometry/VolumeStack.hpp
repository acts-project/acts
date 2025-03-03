// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <vector>

namespace Acts {

/// @class VolumeStack
/// @brief A stack of volumes
/// @note This is a base class for the different types of volume stacks
class VolumeStack : public Volume {
 public:
  /// @param volumes is the vector of volumes
  /// @param direction is the direction of the stack
  /// @param resizeStrategy is the resize strategy
  VolumeStack(std::vector<Volume*>& volumes, AxisDirection direction,
              VolumeResizeStrategy resizeStrategy);

  /// Access the gap volume that were created during attachment or resizing.
  /// @return the vector of gap volumes
  const std::vector<std::shared_ptr<Volume>>& gaps() const;

  /// Check if a volume is a gap volume
  /// @param volume is the volume to check
  /// @return true if the volume is a gap volume, false otherwise
  bool isGapVolume(const Volume& volume) const;

 private:
  /// Helper to get the first volume in the input, and throw an exception if
  /// there is not one.
  /// @param volumes is the vector of volumes
  /// @return the first volume
  static Volume& initialVolume(const std::vector<Volume*>& volumes);

 protected:
  /// @param transform is the transform of the gap volume
  /// @param bounds is the bounds of the gap volume
  /// @return the shared pointer to the gap volume
  std::shared_ptr<Volume> addGapVolume(
      const Transform3& transform, const std::shared_ptr<VolumeBounds>& bounds);

  AxisDirection m_direction{};
  VolumeResizeStrategy m_resizeStrategy{};

  std::vector<std::shared_ptr<Volume>> m_gaps{};
  std::vector<Volume*>& m_volumes;
};
}  // namespace Acts
