// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <iostream>
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @struct ProtoLayer
///
/// Encapsulates min/max boundaries that will be turned into a layer.
/// The struct allows this information to be obtained in a consistent
/// way, or be caller provided.

struct ProtoLayer
{
public:
  double maxX;
  double minX;

  double maxY;
  double minY;

  double maxZ;
  double minZ;

  double maxR;
  double minR;

  double maxPhi;
  double minPhi;

  std::pair<double, double> envX   = {0, 0};
  std::pair<double, double> envY   = {0, 0};
  std::pair<double, double> envZ   = {0, 0};
  std::pair<double, double> envR   = {0, 0};
  std::pair<double, double> envPhi = {0, 0};

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  /// @param surfaces The vector of surfaces to consider
  ProtoLayer(const std::vector<const Surface*>& surfaces);

  /// Constructor
  ///
  /// Loops over a provided vector of surface and calculates the various
  /// min/max values in one go. Also takes into account the thickness
  /// of an associated DetectorElement, if it exists.
  /// @param surfaces The vector of surfaces to consider
  ProtoLayer(const std::vector<std::shared_ptr<const Surface>>& surfaces);

  // normal empty constructor
  ProtoLayer() = default;

  std::ostream&
  dump(std::ostream& sl) const;

  /// Calculates the closest radial distance of a line
  ///
  /// @param pos1 is the first position on the line
  /// @param pos2 is the second position on the line
  ///
  /// @return is the closest distance
  double
  radialDistance(const Vector3D& pos1, const Vector3D& pos2) const;

private:
  /// Helper method which performs the actual min/max calculation
  /// @param surfaces The surfaces to build this protolayer out of
  void
  measure(const std::vector<const Surface*>& surfaces);
};
}  // namespace Acts