// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolumeArrayCreator.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <algorithm>
#include "Acts/Tools/ITrackingVolumeArrayCreator.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class Layer;
class TrackingVolume;

using TrackingVolumeOrderPosition = std::pair<TrackingVolumePtr, Vector3D>;

///@class TrackingVolumeArrayCreator
///
/// The TrackingVolumeArrayCreator is a simple Tool that helps to construct
/// binned arrays of TrackingVolumes for both, confinement in another volume
/// and navigation issues.
///
class TrackingVolumeArrayCreator : public ITrackingVolumeArrayCreator
{
public:
  /// Constructor
  ///
  /// @param logger logging instance
  TrackingVolumeArrayCreator(std::unique_ptr<const Logger> logger
                             = getDefaultLogger("LayerArrayCreator",
                                                Logging::INFO))
    : m_logger(std::move(logger))
  {
  }

  /// Destructor
  ~TrackingVolumeArrayCreator() override = default;

  /// create a tracking volume array
  ///
  /// @param tVolumes is the vector of TrackingVolumes to be
  /// @param bValue is the binning value
  ///
  /// @return new created volume array
  std::shared_ptr<const TrackingVolumeArray>
  trackingVolumeArray(const TrackingVolumeVector& tVolumes,
                      BinningValue                bValue) const override;

  /// Set logging instance
  ///
  /// @param logger is the logging instance to be set
  void
  setLogger(std::unique_ptr<const Logger> logger)
  {
    m_logger = std::move(logger);
  }

private:
  // Private access to the logger method
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<const Logger> m_logger;
};
}