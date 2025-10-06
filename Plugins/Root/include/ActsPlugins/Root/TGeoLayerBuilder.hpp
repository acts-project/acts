// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ILayerBuilder.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/SurfaceBinningMatcher.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/ITGeoIdentifierProvider.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class TGeoMatrix;
class TGeoVolume;
class TGeoNode;

namespace Acts {
class Surface;
class ISurfaceMaterial;
class ITGeoIdentifierProvider;
class LayerCreator;
class ProtoLayerHelper;
}  // namespace Acts

namespace ActsPlugins {

class TGeoDetectorElement;
class ITGeoDetectorElementSplitter;

/// @class TGeoLayerBuilder
///
/// This parses the gGeoManager and looks for a defined combination
/// of volume with contained sensitive detector element. The association
/// is done by matching the names of the TGeoNode / TGeoVolume to the search
/// string.
///
/// The parsing can be restricted to a given parse volume (in r and z),
/// and given some splitting parameters the surfaces can be automatically
/// split into layers.
class TGeoLayerBuilder : public Acts::ILayerBuilder {
 public:
  ///  Helper config structs for volume parsing
  struct LayerConfig {
   public:
    using RangeConfig =
        std::pair<Acts::AxisDirection, std::pair<double, double>>;
    /// Type alias for range configuration with axis direction and bounds

    using SplitConfig = std::pair<Acts::AxisDirection, double>;
    /// Type alias for split configuration with axis direction and position

    /// Identify the search volume by name
    std::string volumeName = "";
    /// Identify the sensor(s) by name
    std::vector<std::string> sensorNames = {};
    /// The local axis definition of TGeo object to Acts::Surface
    std::string localAxes = "XYZ";
    /// Parse ranges: parameter and ranges
    std::vector<RangeConfig> parseRanges = {};
    /// Layer splitting: parameter and tolerance
    std::vector<SplitConfig> splitConfigs = {};
    /// The envelope to be built around the layer
    std::pair<double, double> envelope = {1 * Acts::UnitConstants::mm,
                                          1 * Acts::UnitConstants::mm};
    /// Binning setup in l0: nbins (-1 -> automated), axis binning type
    std::vector<std::pair<int, Acts::BinningType>> binning0 = {
        {-1, Acts::equidistant}};
    /// Binning setup in l1: nbins (-1 -> automated), axis binning type
    std::vector<std::pair<int, Acts::BinningType>> binning1 = {
        {-1, Acts::equidistant}};

    // Default constructor
    LayerConfig() = default;
  };

  /// Type alias for detector element factory function
  using ElementFactory = std::function<std::shared_ptr<TGeoDetectorElement>(
      const TGeoDetectorElement::Identifier&, const TGeoNode&,
      const TGeoMatrix& tGeoMatrix, const std::string& axes, double scalor,
      std::shared_ptr<const Acts::ISurfaceMaterial> material)>;

  /// Default factory function for creating TGeoDetectorElement instances.
  /// @param identifier Detector element identifier
  /// @param tGeoNode ROOT geometry node
  /// @param tGeoMatrix ROOT geometry transformation matrix
  /// @param axes String describing the axes configuration
  /// @param scalor Scaling factor for geometry transformations
  /// @param material Surface material to be assigned to the element
  /// @return Shared pointer to the created TGeoDetectorElement
  static std::shared_ptr<TGeoDetectorElement> defaultElementFactory(
      const TGeoDetectorElement::Identifier& identifier,
      const TGeoNode& tGeoNode, const TGeoMatrix& tGeoMatrix,
      const std::string& axes, double scalor,
      std::shared_ptr<const Acts::ISurfaceMaterial> material);

  /// @struct Config
  /// @brief nested configuration struct for steering of the layer builder
  struct Config {
    /// String based identification
    std::string configurationName = "undefined";
    /// Unit conversion
    double unit = 1 * Acts::UnitConstants::cm;
    /// Create an identifier from TGeoNode
    std::shared_ptr<const ITGeoIdentifierProvider> identifierProvider = nullptr;
    /// Split TGeoElement if a splitter is provided
    std::shared_ptr<const ITGeoDetectorElementSplitter>
        detectorElementSplitter = nullptr;
    /// Factory for creating detector elements based on TGeoNodes
    ElementFactory detectorElementFactory = defaultElementFactory;
    /// Layer creator
    std::shared_ptr<const Acts::LayerCreator> layerCreator = nullptr;
    /// ProtoLayer helper
    std::shared_ptr<const Acts::ProtoLayerHelper> protoLayerHelper = nullptr;
    /// Configuration is always | n | c | p |
    std::array<std::vector<LayerConfig>, 3> layerConfigurations;
    /// Split tolerances in R
    std::array<double, 3> layerSplitToleranceR = {-1., -1., -1.};
    /// Split tolerances in Z
    std::array<double, 3> layerSplitToleranceZ = {-1., -1., -1.};
    /// Automated binning & tolerances
    bool autoSurfaceBinning = false;
    /// The surface binning matcher
    Acts::SurfaceBinningMatcher surfaceBinMatcher;
  };

  /// Constructor
  /// @param config is the configuration struct
  /// @param logger the local logging instance
  explicit TGeoLayerBuilder(const Config& config,
                            std::unique_ptr<const Acts::Logger> logger =
                                Acts::getDefaultLogger("TGeoLayerBuilder",
                                                       Acts::Logging::INFO));

  /// Destructor
  ~TGeoLayerBuilder() override;

  /// LayerBuilder interface method - returning the layers at negative side
  ///
  /// @param gctx the geometry context for this build call
  ///
  /// @return Vector of negative-side layers
  const Acts::LayerVector negativeLayers(
      const Acts::GeometryContext& gctx) const final;

  /// LayerBuilder interface method - returning the central layers
  ///
  /// @param gctx the geometry context for this build call
  ///
  /// @return Vector of central layers
  const Acts::LayerVector centralLayers(
      const Acts::GeometryContext& gctx) const final;

  /// LayerBuilder interface method - returning the layers at positive side
  ///
  /// @param gctx the geometry context for this build call
  ///
  /// @return Vector of positive-side layers
  const Acts::LayerVector positiveLayers(
      const Acts::GeometryContext& gctx) const final;

  /// Name identification
  /// @return Reference to the configuration name string
  const std::string& identification() const final;

  /// Set the configuration object
  /// @param config is the configuration struct
  void setConfiguration(const Config& config);

  /// Get the configuration object
  /// @return Copy of the current configuration
  Config getConfiguration() const;

  /// Set logging instance
  /// @param newLogger Logger instance to use for debug output
  void setLogger(std::unique_ptr<const Acts::Logger> newLogger);

  /// Return the created detector elements
  /// @return Reference to vector of detector elements
  const std::vector<std::shared_ptr<const TGeoDetectorElement>>&
  detectorElements() const;

 private:
  /// Configuration object
  Config m_cfg;

  /// layer types
  std::array<std::string, 3> m_layerTypes = {"Negative", "Central", "Positive"};

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// @todo make clear where the TGeoDetectorElement lives
  std::vector<std::shared_ptr<const TGeoDetectorElement>> m_elementStore;

  /// Private helper method : build layers
  ///
  /// @param gctx the geometry context of this call
  /// @param layers is goint to be filled
  /// @param type is the indication which ones to build -1 | 0 | 1
  void buildLayers(const Acts::GeometryContext& gctx, Acts::LayerVector& layers,
                   int type = 0);

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
  for (const auto& splitPar : parameters) {
    if (std::abs(test - splitPar) < tolerance) {
      found = true;
    }
  }
  if (!found) {
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

}  // namespace ActsPlugins
