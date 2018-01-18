// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_LAYERS_PROTOLAYER_H
#define ACTS_LAYERS_PROTOLAYER_H 1

#include <iostream>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"

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
  ProtoLayer(std::vector<const Surface*> surfaces);

  // normal empty constructor
  ProtoLayer(){};

  std::ostream&
  dump(std::ostream& sl) const;

private:
  /// Calculates the closest radial distance of a line
  ///
  /// @param pos1 is the first position on the line
  /// @param pos2 is the second position on the line
  ///
  /// @return is the closest distance
  double
  radialDistance(const Vector3D& pos1, const Vector3D& pos2);
};
}  // namespace Acts

#endif  // ACTS_LAYERS_PROTOLAYER_H
