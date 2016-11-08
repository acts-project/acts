// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingVolumeBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TOOLS_ITRACKINGVOLUMEBUILDER_H
#define ACTS_TOOLS_ITRACKINGVOLUMEBUILDER_H 1

#include <memory>
#include <tuple>
#include <vector>

namespace Acts {

class VolumeBounds;
class TrackingVolume;
class Layer;
class Volume;
typedef std::shared_ptr<const TrackingVolume> TrackingVolumePtr;
typedef std::shared_ptr<const VolumeBounds>   VolumeBoundsPtr;
typedef std::shared_ptr<const Layer>          LayerPtr;
typedef std::vector<LayerPtr>                 LayerVector;
///  @class ITrackingVolumeBuilder
///
/// Interface class ITrackingVolumeBuilders
///
/// this returns the sub-detector tracking volume that is wrapped by the next
/// outer one
/// in the TrackingGeometry building process
///
/// If an innerVolume is given, this is wrapped
/// If a VolumeBounds object is given this defines the maximum extent.
///

class ITrackingVolumeBuilder
{
public:
  /// Virtual destructor
  virtual ~ITrackingVolumeBuilder() = default;

  /// ITrackingVolumeBuilder interface method
  ///
  /// @param insideVolume is an (optional) volume to be wrapped
  /// @param outsideBounds is an (optional) outside confinement
  ///
  /// @return shared pointer to a newly created TrackingVolume
  virtual TrackingVolumePtr
  trackingVolume(TrackingVolumePtr insideVolume  = nullptr,
                 VolumeBoundsPtr   outsideBounds = nullptr) const = 0;
};

}  // end of namespace

#endif  // ACTS_TOOLS_ITRACKINGVOLUMEBUILDER_H
