// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoLayerBuilder.h, ACTS project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_TGEOPLUGINS_TGEOLAYERBUILDER_H
#define ACTS_TGEOPLUGINS_TGEOLAYERBUILDER_H

#include "ACTS/Tools/ILayerBuilder.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Logger.hpp"

class TGeoMatrix;
class TGeoVolume;
class TGeoNode;

namespace Acts {

class ILayerCreator;
class TGeoDetectorElement;
class Surface;

typedef std::pair<TGeoNode*, std::shared_ptr<Transform3D>> NodeTransform;

/// @class TGeoLayerBuilder
/// works on the gGeoManager, as this is filled from GDML
class TGeoLayerBuilder : public ILayerBuilder
{
public:
  ///  Helper config structs for volume parsin
  struct LayerConfig
  {
  public:
    /// identify the layer by name
    std::string layerName;
    /// identify the sensor by name
    std::string sensorName;
    // the local axis definition
    std::string localAxes;
    // the envolpoe
    std::pair<double, double> envelope;
    /// define the number of bins in loc0
    size_t binsLoc0;
    /// define the number of bins in loc1
    size_t binsLoc1;

    LayerConfig()
      : layerName("")
      , sensorName("")
      , localAxes("xyz")
      , envelope(std::pair<double, double>(0., 0.))
      , binsLoc0(0)
      , binsLoc1(0)
    {
    }
  };

  /// @struct Config
  /// nested configuration struct for steering of the layer builder
  struct Config
  {
    /// string based identification
    std::string configurationName = "undefined";
    // unit conversion
    double unit = 10;
    // set visibility flag
    bool setVisibility;
    // layer creator
    std::shared_ptr<ILayerCreator> layerCreator = nullptr;
    // configurations
    std::vector<LayerConfig> negativeLayerConfigs;
    std::vector<LayerConfig> centralLayerConfigs;
    std::vector<LayerConfig> positiveLayerConfigs;
  };

  /// Constructor
  /// @param cfg is the configuration struct
  /// @param logger the local logging instance
  TGeoLayerBuilder(const Config&           cfg,
                   std::unique_ptr<Logger> logger
                   = getDefaultLogger("LayerArrayCreator", Logging::INFO));

  /// Destructor
  ~TGeoLayerBuilder();

  /// LayerBuilder interface method - returning the layers at negative side
  const LayerVector
  negativeLayers() const final;

  /// LayerBuilder interface method - returning the central layers
  const LayerVector
  centralLayers() const final;

  /// LayerBuilder interface method - returning the layers at negative side
  const LayerVector
  positiveLayers() const final;

  /// Name identification
  const std::string&
  identification() const final;

  /// set the configuration object
  /// @param cfg is the configuration struct
  void
  setConfiguration(const Config& cfg);

  /// get the configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<Logger> logger);

private:
  /// configruation object
  Config m_cfg;

  /// Private access to the logger
  const Logger&
  logger() const
  {
    return *m_logger;
  }

  /// logging instance
  std::unique_ptr<Logger> m_logger;

  /// @todo make clear where the TGeoDetectorElement lives
  mutable std::vector<std::shared_ptr<TGeoDetectorElement>> m_elementStore;

  /// Private helper function to parse the geometry tree
  void
  collectSensitive(std::vector<const Surface*>& layerSurfaces,
                   TGeoVolume*                  tgVolume,
                   TGeoNode*                    tgNode,
                   const TGeoMatrix&            ctGlobal,
                   const LayerConfig&           layerConfig,
                   int                          type,
                   bool                         correctVolume = false,
                   const std::string&           offset        = "") const;

  // Private helper mehtod : build layers
  // @param layers is goint to be filled
  // @param type is the indication which ones to build -1 | 0 | 1
  void
  buildLayers(LayerVector& layers, int type = 0) const;
};

inline TGeoLayerBuilder::Config
TGeoLayerBuilder::getConfiguration() const
{
  return m_cfg;
}

inline const std::string&
TGeoLayerBuilder::identification() const
{
  return m_cfg.configurationName;
}
}

#endif
