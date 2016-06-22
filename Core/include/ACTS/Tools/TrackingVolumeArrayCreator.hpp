// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolumeArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_TRACKINGVOLUMEARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_TRACKINGVOLUMEARRAYCREATOR_H

#include "ACTS/Tools/ITrackingVolumeArrayCreator.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include <algorithm>
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

class Layer;
class TrackingVolume;

typedef std::pair<TrackingVolumePtr, Vector3D> TrackingVolumeOrderPosition;

///@class TrackingVolumeArrayCreator
///
/// The TrackingVolumeArrayCreator is a simple Tool that helps to construct
/// binned arrays of TrackingVolumes for both, confinement in another volume
/// and navigation issues.
///
class TrackingVolumeArrayCreator : public ITrackingVolumeArrayCreator
{
public:
  /// @struct Config
  /// Nested configuration struct for this array creator
  struct Config
  {
    std::shared_ptr<Logger> logger;  ///< logging instance

    Config() : logger(getDefaultLogger("LayerArrayCreator", Logging::INFO)) {}
  };

  /// Constructor 
  /// @param cfg is the configuration class
  TrackingVolumeArrayCreator(const Config& cfg) : m_cfg(cfg) {}
  
  /// Destructor 
  virtual ~TrackingVolumeArrayCreator() = default;

  /// create a tracking volume array 
  /// @param vols is the vector of TrackingVolumes to be 
  /// @param bVals is the binning value
  std::shared_ptr<const TrackingVolumeArray>
  trackingVolumeArray(const TrackingVolumeVector& vols,
                      BinningValue                bVal) const;

  /// Set the configuration
  /// @param cfg configuration struct to be set
  void
  setConfiguration(const Config& cfg)
  {
    m_cfg = cfg;
  }

  /// Get configuration method 
  Config
  getConfiguration() const
  {
    return m_cfg;
  }

private:
  /// Configuration struct 
  Config m_cfg;
  
  // Private access to the logger method
  const Logger&
  logger() const
  {
    return *m_cfg.logger;
  }
};
}

#endif
