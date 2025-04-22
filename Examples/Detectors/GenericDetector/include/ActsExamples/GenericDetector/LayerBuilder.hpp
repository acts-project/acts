// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ILayerBuilder.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"
#include "ActsExamples/GenericDetector/ProtoLayerCreator.hpp"

namespace ActsExamples::Generic {

/// @class LayerBuilder
///
/// The LayerBuilder is able to build cylinder & disc layers
/// from hard-coded input.
/// This is meant for the simple detector examples.
class LayerBuilder : public Acts::ILayerBuilder {
 public:
  /// @struct Config
  /// Nested configuration struct for the LayerBuilderT
  struct Config {
    /// The string based identification
    std::string layerIdentification = "";

    /// The pre-produced proto layers for the central part
    std::vector<ProtoLayerSurfaces> centralProtoLayers;

    /// the material concentration: -1 inner, 0 central, 1 outer
    std::vector<int> centralLayerMaterialConcentration;
    /// The assigned surface material
    std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>
        centralLayerMaterial;

    /// The pre-produced proto layers for the negative part
    std::vector<ProtoLayerSurfaces> negativeProtoLayers;

    /// The pre-produced proto layers for the positive part
    std::vector<ProtoLayerSurfaces> positiveProtoLayers;

    /// The material concentration: -1 inner, 0 central, 1 outer
    std::vector<int> posnegLayerMaterialConcentration;

    /// The surface material
    std::vector<std::shared_ptr<const Acts::ISurfaceMaterial>>
        posnegLayerMaterial;

    /// Helper tools: layer creator
    std::shared_ptr<const Acts::LayerCreator> layerCreator = nullptr;
    /// Helper tools: central passive layer builder
    std::shared_ptr<const Acts::ILayerBuilder> centralPassiveLayerBuilder =
        nullptr;
    /// Helper tools: p/n passive layer builder
    std::shared_ptr<const Acts::ILayerBuilder> posnegPassiveLayerBuilder =
        nullptr;
  };

  /// Constructor
  /// @param glbConfig is the configuration class
  explicit LayerBuilder(const Config& cfg,
                        std::unique_ptr<const Acts::Logger> logger =
                            Acts::getDefaultLogger("LayerBuilder",
                                                   Acts::Logging::INFO));

  /// LayerBuilder interface method - returning the layers at negative side
  const Acts::LayerVector negativeLayers(
      const Acts::GeometryContext& gctx) const override;

  /// LayerBuilder interface method - returning the central layers
  const Acts::LayerVector centralLayers(
      const Acts::GeometryContext& gctx) const override;

  /// LayerBuilder interface method - returning the layers at positive side
  const Acts::LayerVector positiveLayers(
      const Acts::GeometryContext& gctx) const override;

  /// ILayerBuilder method
  const std::string& identification() const override {
    return m_cfg.layerIdentification;
  }

 private:
  const Acts::LayerVector constructEndcapLayers(
      const Acts::GeometryContext& gctx, int side) const;

  /// Configuration member
  Config m_cfg;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// the logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples::Generic
