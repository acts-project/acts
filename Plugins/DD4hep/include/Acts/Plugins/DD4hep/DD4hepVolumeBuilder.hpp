// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tools/IConfinedTrackingVolumeBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "Acts/Detector/TrackingVolume.hpp"

class TrackingVolume;
using MutableTrackingVolumePtr    = std::shared_ptr<TrackingVolume>;
using MutableTrackingVolumeVector = std::vector<MutableTrackingVolumePtr>;

class TGeoMatrix;

namespace dd4hep {
class DetElement;
}

namespace Acts {

/// @brief build layers of one cylinder-endcap setup from DD4hep input
///
/// This class is an implementation of the Acts::ILayerBuilder,
/// creating the central (layers of barrel), the negative and positive layers
/// (layers of endcaps) of one hierarchy (e.g. PixelDetector, StripDetector,...)
/// with input from DD4hep.

class DD4hepVolumeBuilder : public IConfinedTrackingVolumeBuilder
{
public:
  /// @struct Config
  /// nested configuration struct for steering of the layer builder
  struct Config
  {
    /// string based identification
    std::string configurationName = "undefined";

    std::vector<dd4hep::DetElement> centralVolumes;
  };

  /// Constructor
  /// @param config is the configuration struct
  /// @param logger is the logging instance
  DD4hepVolumeBuilder(const Acts::DD4hepVolumeBuilder::Config& config,
                      std::unique_ptr<const Logger>            logger);
  /// Destructor
  ~DD4hepVolumeBuilder() override;

  //~ /// LayerBuilder interface method
  //~ /// @return  the layers at negative side
  //~ const LayerVector
  //~ negativeLayers() const final;

  //~ /// LayerBuilder interface method
  //~ /// @return the layers at the central sector
  MutableTrackingVolumeVector
  centralVolumes() const final;

  //~ /// LayerBuilder interface method
  //~ /// @return  the layers at positive side
  //~ const LayerVector
  //~ positiveLayers() const final;

  /// Name identification
  /// @return the string based identification of this configuration
  const std::string&
  identification() const final;

  /// set the configuration object
  /// @param config is the configuration struct
  void
  setConfiguration(const Config& config);

  /// get the configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<const Logger> logger);

private:
  /// configruation object
  Config m_cfg;

  /// logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  std::shared_ptr<const Acts::Transform3D>
  convertTransform(const TGeoMatrix* tGeoTrans) const;
};

inline const std::string&
DD4hepVolumeBuilder::identification() const
{
  return m_cfg.configurationName;
}

inline DD4hepVolumeBuilder::Config
DD4hepVolumeBuilder::getConfiguration() const
{
  return m_cfg;
}

}  // namespace
