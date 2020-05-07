// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TGeoLayerBuilder.h, Acts project, TGeoDetector plugin
///////////////////////////////////////////////////////////////////

#pragma once
#include <climits>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ILayerBuilder.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Plugins/TGeo/ITGeoIdentifierProvider.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

class TGeoMatrix;
class TGeoVolume;
class TGeoNode;

namespace Acts {

class TGeoDetectorElement;
class Surface;

using NodeTransform = std::pair<TGeoNode*, std::shared_ptr<const Transform3D>>;
using namespace Acts::UnitLiterals;

/// @class TGeoLayerBuilder
///
/// This parses the gGeoManager and looks for a defined combination
/// of volume with contained sensitive detector element. The association
/// is done by matching the names of the TGeoNode / TGeoVolume to the search
/// string.
///
/// The parsing can be restricted to a given parse volume (in r and z),
/// and given some splitting parameters the surfaces can be automatically be
/// split into layers.
class TGeoLayerBuilder : public ILayerBuilder {
 public:
  ///  Helper config structs for volume parsin
  struct LayerConfig {
   public:
    /// Identify the layer by name
    std::string layerName = "";
    /// Identify the sensor by name
    std::string sensorName = "";
    /// The local axis definition of TGeo object to Acts::Surface
    std::string localAxes = "xyz";
    /// Parse area :  in r
    std::pair<double, double> parseRangeR = {
        0., std::numeric_limits<double>::max()};
    /// Parse area : in z
    std::pair<double, double> parseRangeZ = {
        -std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()};
    /// The envelope to be built around the layer
    std::pair<double, double> envelope = {0_mm, 0_mm};
    /// Layer splitting: r min/max of the split range
    std::pair<double, double> splitRangeR = {
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max()};
    /// Layer splitting: r the result of the splitting
    std::vector<double> splitParametersR = {};
    /// Layer splitting: z min/max of the split range
    std::pair<double, double> splitRangeZ = {
        std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max()};
    /// Layer splitting: z the result of the splitting
    std::vector<double> splitParametersZ = {};
    /// Define the number of bins in loc0
    size_t binsLoc0{1};
    /// Define the number of bins in loc1
    size_t binsLoc1{1};

    // Default constructor
    LayerConfig()
        : layerName(""),
          sensorName(""),
          localAxes("XZY"),
          envelope(std::pair<double, double>(1_mm, 1_mm)) {}
  };

  /// @struct Config
  /// @brief nested configuration struct for steering of the layer builder
  struct Config {
    /// String based identification
    std::string configurationName = "undefined";
    /// Unit conversion
    double unit = 1_cm;
    /// Create an indentifier from TGeoNode
    std::shared_ptr<const ITGeoIdentifierProvider> identifierProvider = nullptr;
    /// Layer creator
    std::shared_ptr<const LayerCreator> layerCreator = nullptr;
    /// Configuration is always | n | c | p |
    std::array<std::vector<LayerConfig>, 3> layerConfigurations;
    /// Split tolerances in R
    std::array<double, 3> layerSplitToleranceR = {-1., -1., -1.};
    /// Split tolerances in Z
    std::array<double, 3> layerSplitToleranceZ = {-1., -1., -1.};
    /// Check for ring layout when building a volume
    bool checkRingLayout = false;
    /// Tolerance for ring detection and association
    double ringTolerance = 0_mm;
    /// Special debug output, is very verbose and hence needs
    /// an additional switch to log level
    bool nodeSearchDebug = false;
  };

  /// Constructor
  /// @param config is the configuration struct
  /// @param logger the local logging instance
  TGeoLayerBuilder(const Config& config,
                   std::unique_ptr<const Logger> logger =
                       getDefaultLogger("LayerArrayCreator", Logging::INFO));

  /// Destructor
  ~TGeoLayerBuilder() override;

  /// LayerBuilder interface method - returning the layers at negative side
  ///
  /// @param gctx the geometry context for this build call
  ///
  const LayerVector negativeLayers(const GeometryContext& gctx) const final;

  /// LayerBuilder interface method - returning the central layers
  ///
  /// @param gctx the geometry context for this build call
  ///
  const LayerVector centralLayers(const GeometryContext& gctx) const final;

  /// LayerBuilder interface method - returning the layers at negative side
  ///
  /// @param gctx the geometry context for this build call
  ///
  const LayerVector positiveLayers(const GeometryContext& gctx) const final;

  /// Name identification
  const std::string& identification() const final;

  /// Set the configuration object
  /// @param config is the configuration struct
  void setConfiguration(const Config& config);

  /// Get the configuration object
  Config getConfiguration() const;

  /// Set logging instance
  void setLogger(std::unique_ptr<const Logger> newLogger);

  /// Return the created detector elements
  const std::vector<std::shared_ptr<const TGeoDetectorElement>>&
  detectorElements() const;

 private:
  /// Configuration object
  Config m_cfg;

  /// layer types
  std::array<std::string, 3> m_layerTypes = {"Negative", "Central", "Positive"};

  /// Private access to the logger
  const Logger& logger() const { return *m_logger; }

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// @todo make clear where the TGeoDetectorElement lives
  std::vector<std::shared_ptr<const TGeoDetectorElement>> m_elementStore;

  /// Private helper function to parse the geometry tree
  ///
  /// @param gcts the geometry context of this call
  /// @param layerSurfaces are the surfaces that build the layer
  /// @param tgVolume is the current volume
  /// @param tgNode the current Node (branch)
  /// @param tgTransform is the current relative transform
  /// @param layerConfig is the configuration to be filled
  /// @param type is ( n | c | p ) as of  ( -1 | 0 | 1 )
  /// @param correctBranch is the branch hit
  /// @param offset is a string offset for the screen output
  void resolveSensitive(
      const GeometryContext& gctx,
      std::vector<std::shared_ptr<const Surface>>& layerSurfaces,
      TGeoVolume* tgVolume, TGeoNode* tgNode, const TGeoMatrix& tgTransform,
      LayerConfig& layerConfig, int type, bool correctBranch = false,
      const std::string& offset = "");

  /// Private helper method : build layers
  ///
  /// @param gcts the geometry context of this call
  /// @param layers is goint to be filled
  /// @param type is the indication which ones to build -1 | 0 | 1
  void buildLayers(const GeometryContext& gctx, LayerVector& layers,
                   int type = 0);

  /// Private helper method : match string with wildcards
  /// @param wc is the one with the potential wildcard
  /// @param test is the test string
  bool match(const char* first, const char* second) const;

  /// Private helper method : register splitting input
  void registerSplit(std::vector<double>& parameters, double test,
                     double tolerance, std::pair<double, double>& range) const;
};

inline void TGeoLayerBuilder::registerSplit(
    std::vector<double>& parameters, double test, double tolerance,
    std::pair<double, double>& range) const {
  bool found = false;
  // min/max setting
  range.first = std::min(range.first, test);
  range.second = std::max(range.second, test);
  // Loop and find the split parameters
  for (auto& splitPar : parameters) {
    if (std::abs(test - splitPar) < tolerance) {
      found = true;
    }
  }
  if (not found) {
    parameters.push_back(test);
  }
}

inline TGeoLayerBuilder::Config TGeoLayerBuilder::getConfiguration() const {
  return m_cfg;
}

inline const std::vector<std::shared_ptr<const TGeoDetectorElement>>&
TGeoLayerBuilder::detectorElements() const {
  return m_elementStore;
}

inline const std::string& TGeoLayerBuilder::identification() const {
  return m_cfg.configurationName;
}

// The main function that checks if two given strings
// match. The first string may contain wildcard characters
inline bool TGeoLayerBuilder::match(const char* first,
                                    const char* second) const {
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
}  // namespace Acts
