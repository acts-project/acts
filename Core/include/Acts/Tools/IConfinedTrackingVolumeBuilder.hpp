// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
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
using MutableTrackingVolumePtr    = std::shared_ptr<TrackingVolume>;
using MutableTrackingVolumeVector = std::vector<MutableTrackingVolumePtr>;

class IConfinedTrackingVolumeBuilder
{
public:
  /// Virtual destructor
  virtual ~IConfinedTrackingVolumeBuilder() = default;

  virtual MutableTrackingVolumeVector
  centralVolumes() const = 0;

  virtual const std::string&
  identification() const = 0;
};

}  // namespace
