// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingVolumeBuilder.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <memory>
#include <vector>

namespace Acts {

class TrackingVolume;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;
using MutableTrackingVolumeVector = std::vector<MutableTrackingVolumePtr>;

/// @brief This is an interface class for constructing TrackingVolumes whose are
/// confined in a mother-TrackingVolume
class IConfinedTrackingVolumeBuilder {
 public:
  /// Virtual destructor
  virtual ~IConfinedTrackingVolumeBuilder() = default;

  /// Interface for constructing a vector of confined TrackingVolumes
  virtual MutableTrackingVolumeVector centralVolumes() const = 0;

  /// Interface for retreiving the identification string of the confined volumes
  virtual const std::string& identification() const = 0;
};

}  // namespace Acts