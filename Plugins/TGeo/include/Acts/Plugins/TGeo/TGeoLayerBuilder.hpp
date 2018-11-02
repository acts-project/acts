// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoLayerBuilder.h, Acts project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Tools/ILayerBuilder.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

class TGeoMatrix;
class TGeoVolume;
class TGeoNode;

namespace Acts {

class TGeoDetectorElement;
class Surface;

using NodeTransform = std::pair<TGeoNode*, std::shared_ptr<const Transform3D>>;

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
    std::string layerName = "";
    /// identify the sensor by name
    std::string sensorName = "";
    // the local axis definition
    std::string localAxes = "xyz";
    // the envolpoe
    std::pair<double, double> envelope;
    /// define the number of bins in loc0
    size_t binsLoc0{100};
    /// define the number of bins in loc1
    size_t binsLoc1{100};

    LayerConfig()
      : layerName("")
      , sensorName("")
      , localAxes("XZY")
      , envelope(std::pair<double, double>(1., 1.))
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
    std::shared_ptr<const LayerCreator> layerCreator = nullptr;
    // configurations
    std::vector<LayerConfig> negativeLayerConfigs;
    std::vector<LayerConfig> centralLayerConfigs;
    std::vector<LayerConfig> positiveLayerConfigs;
  };

  /// Constructor
  /// @param config is the configuration struct
  /// @param logger the local logging instance
  TGeoLayerBuilder(const Config&                 config,
                   std::unique_ptr<const Logger> logger
                   = getDefaultLogger("LayerArrayCreator", Logging::INFO));

  /// Destructor
  ~TGeoLayerBuilder() override;

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
  /// @param config is the configuration struct
  void
  setConfiguration(const Config& config);

  /// get the configuration object
  Config
  getConfiguration() const;

  /// set logging instance
  void
  setLogger(std::unique_ptr<const Logger> newLogger);

  /// Return the created detector elements
  const std::vector<std::shared_ptr<const TGeoDetectorElement>>&
  detectorElements() const;

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
  std::unique_ptr<const Logger> m_logger;

  /// @todo make clear where the TGeoDetectorElement lives
  std::vector<std::shared_ptr<const TGeoDetectorElement>> m_elementStore;

  /// Private helper function to parse the geometry tree
  void
  resolveSensitive(std::vector<std::shared_ptr<const Surface>>& layerSurfaces,
                   TGeoVolume*                                  tgVolume,
                   TGeoNode*                                    tgNode,
                   const TGeoMatrix&                            tgTransform,
                   const LayerConfig&                           layerConfig,
                   int                                          type,
                   bool               correctBranch = false,
                   const std::string& offset        = "");

  // Private helper method : build layers
  // @param layers is goint to be filled
  // @param type is the indication which ones to build -1 | 0 | 1
  void
  buildLayers(LayerVector& layers, int type = 0);

  // Private helper method : match string with wildcards
  // @param wc is the one with the potential wildcard
  // @param test is the test string
  bool
  match(const char* first, const char* second) const;
};

inline TGeoLayerBuilder::Config
TGeoLayerBuilder::getConfiguration() const
{
  return m_cfg;
}

inline const std::vector<std::shared_ptr<const TGeoDetectorElement>>&
TGeoLayerBuilder::detectorElements() const
{
  return m_elementStore;
}

inline const std::string&
TGeoLayerBuilder::identification() const
{
  return m_cfg.configurationName;
}

// The main function that checks if two given strings
// match. The first string may contain wildcard characters
inline bool
TGeoLayerBuilder::match(const char* first, const char* second) const
{
  // If we reach at the end of both strings, we are done
  if (*first == '\0' && *second == '\0') {
    return true;
  }

  // Make sure that the characters after '*' are present
  // in second string. This function assumes that the first
  // string will not contain two consecutive '*'
  if (*first == '*' && *(first + 1) != '\0' && *second == '\0') {
    return false;
  }

  // If the first string contains '?', or current characters
  // of both strings match
  if (*first == '?' || *first == *second) {
    return match(first + 1, second + 1);
  }

  // If there is *, then there are two possibilities
  // a) We consider current character of second string
  // b) We ignore current character of second string.
  if (*first == '*') {
    return match(first + 1, second) || match(first, second + 1);
  }
  return false;
}
}
