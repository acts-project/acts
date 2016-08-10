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
#include "ACTS/Utilities/Logger.hpp"

class TGeoVolume;
class TGeoNode;

namespace Acts {

class ILayerCreator;
class TGeoDetectorElement;

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
    /// define the number of bins in loc0
    size_t binsLoc0;
    /// define the number of bins in loc1
    size_t binsLoc1;

    LayerConfig() : layerName(""), sensorName(""), binsLoc0(0), binsLoc1(0) {}
  };

  /// @struct Config
  /// nested configuration struct for steering of the layer builder
  class Config
  {
  public:
    /// logging instance 
    std::shared_ptr<Logger>        logger; 
    /// string based identification
    std::string                    configurationName;
    // layer creator
    std::shared_ptr<ILayerCreator> layerCreator;
    // configurations 
    std::vector<LayerConfig>       negativeLayerConfigs;
    std::vector<LayerConfig>       centralLayerConfigs;
    std::vector<LayerConfig>       positiveLayerConfigs;

    Config()
      : logger(getDefaultLogger("LayerArrayCreator", Logging::INFO))
      , configurationName("Undefined")
      , layerCreator(nullptr)
      , negativeLayerConfigs()
      , centralLayerConfigs()
      , positiveLayerConfigs()
    {
    }
  };

  /// Constructor 
  /// @param cfg is the configuration struct
  TGeoLayerBuilder(const Config& cfg);

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

private:
  /// configruation object
  Config m_cfg;

  /// @TODO make clear where the TGeoDetectorElement lives
  mutable std::vector<std::shared_ptr<TGeoDetectorElement>> m_elementStore;

  /// Private helper function to parse the geometry tree
  void
  collectSensitive(TGeoVolume*             tgVolume,
                   TGeoNode*               tgNode,
                   const std::string&      sensitiveName,
                   std::vector<TGeoNode*>& collectedVolumes) const;
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
