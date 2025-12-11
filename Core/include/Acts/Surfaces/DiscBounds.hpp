// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

/// @class DiscBounds
///
/// @image html DiscBounds.gif
///
/// common base class for all bounds that are in a r/phi frame
///  - simply introduced to avoid wrong bound assignments to surfaces
///
class DiscBounds : public SurfaceBounds {
 public:
  /// Return method for inner Radius
  /// @return Minimum radius of the disc bounds
  virtual double rMin() const = 0;

  /// Return method for outer Radius
  /// @return Maximum radius of the disc bounds
  virtual double rMax() const = 0;

  /// Return the vertices
  ///
  /// @param quarterSegments The number of segments used to describe a quarter
  /// of a circle, if it is 1, then only the extrema points in phi are inserted
  /// next to the segment corners
  ///
  /// @return vector for vertices in 2D
  virtual std::vector<Vector2> vertices(
      unsigned int quarterSegments = 2u) const = 0;

  /// Returns a reference radius for binning
  /// @return Reference radius value used for binning purposes
  virtual double binningValueR() const = 0;

  /// Returns a refererance phi for binning
  /// @return Reference phi angle value used for binning purposes
  virtual double binningValuePhi() const = 0;

  /// Returns true for full phi coverage
  /// @return True if the bounds cover the full azimuthal range, false otherwise
  virtual bool coversFullAzimuth() const = 0;

  /// Checks if it's inside the radius
  /// @param R The radius value to check
  /// @param tolerance The tolerance for the check
  /// @return True if the radius is within the bounds, false otherwise
  virtual bool insideRadialBounds(double R, double tolerance = 0.) const = 0;
};

}  // namespace Acts
