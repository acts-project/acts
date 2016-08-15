// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingVolumeArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEARRAYCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ITRACKINGVOLUMEARRAYCREATOR_H 1

#include <memory>
#include <vector>
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/BinningType.hpp"

namespace Acts {

class TrackingVolume;
typedef std::shared_ptr<const TrackingVolume> TrackingVolumePtr;

/// @typedef TrackingVolumeArray
typedef BinnedArray<TrackingVolumePtr> TrackingVolumeArray;
typedef std::vector<TrackingVolumePtr> TrackingVolumeVector;

/// @class ITrackingVolumeArrayCreator
///
/// Interface class ITrackingVolumeArrayCreators It inherits from IAlgTool.
///
/// It is designed to centralize the code to create
/// Arrays of Tracking Volumes for both:
///
///   - confinement in another TrackingVolume
///   - navigation and glueing
///
/// Arrays for glueing and confinement are often the same,
/// therefore the newly created TrackingVolumeArray is done by a shared_ptr
class ITrackingVolumeArrayCreator
{
public:
  /// Virtual destructor
  virtual ~ITrackingVolumeArrayCreator() = default;

  /// TrackingVolumeArrayCreator interface method - creates array depending on
  /// the binning type
  /// @param vols are the TrackingVolumes ordered in a tracker
  /// @param bVal is the binning value for the volume binning
  virtual std::shared_ptr<const TrackingVolumeArray>
  trackingVolumeArray(const TrackingVolumeVector& vols,
                      BinningValue                bVal) const = 0;
};
}  // end of namespace

#endif
