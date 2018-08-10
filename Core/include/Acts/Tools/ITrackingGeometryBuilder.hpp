// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingGeometryBuilder.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>

namespace Acts {
class TrackingGeometry;

/// @class ITrackingGeometryBuilder
///
/// Interface class for the TrackingGeometry building,
/// this is used by the TrackingGeometrySvc to build the geoemtry.
///
/// The TrackingGeometry is written to the detector store and thus not created
/// as a std::shared_ptr.
///
/// The TrackingGeometry is returned as a non-const object in order to recreate
/// from conditions callback if necessary.
///
class ITrackingGeometryBuilder
{
public:
  /// Virtual destructor
  virtual ~ITrackingGeometryBuilder() = default;

  /// TrackingGeometry Interface methode
  /// @return unique pointer to a newly created TrackingGeometry
  virtual std::unique_ptr<const TrackingGeometry>
  trackingGeometry() const = 0;
};
}  // namespace