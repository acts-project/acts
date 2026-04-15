// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"

namespace Acts {
class GeometryContext;
}

namespace ActsPlugins {
class DD4hepDetectorElement;
}

namespace ActsExamples {

class OpenDataDetector final : public DD4hepDetectorBase {
 public:
  struct Config : public DD4hepDetectorBase::Config {
    enum class ConstructionMethod {
      BarrelEndcap,
      DirectLayer,
      DirectLayerGrouped
    };

    using ElementFactory =
        std::function<std::shared_ptr<ActsPlugins::DD4hepDetectorElement>(
            const dd4hep::DetElement& element, ActsPlugins::TGeoAxes axes,
            double scale)>;

    ElementFactory detectorElementFactory = defaultDetectorElementFactory;

    /// Select the conversion style used to construct the Gen3 ODD geometry.
    ConstructionMethod constructionMethod = ConstructionMethod::BarrelEndcap;

    /// Envelope for the blueprint root (world volume). Values in mm.
    Acts::ExtentEnvelope blueprintEnvelope =
        Acts::ExtentEnvelope::Zero()
            .set(Acts::AxisDirection::AxisZ, {20., 20.})
            .set(Acts::AxisDirection::AxisR, {0., 20.});

    /// Envelope for layer volumes. Values in mm.
    Acts::ExtentEnvelope layerEnvelope =
        Acts::ExtentEnvelope::Zero()
            .set(Acts::AxisDirection::AxisZ, {2., 2.})
            .set(Acts::AxisDirection::AxisR, {2., 2.});
  };

  static std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
  defaultDetectorElementFactory(const dd4hep::DetElement& element,
                                ActsPlugins::TGeoAxes axes, double scale);

  explicit OpenDataDetector(const Config& cfg,
                            const Acts::GeometryContext& gctx);

  const Config& config() const override;

 private:
  void construct(const Acts::GeometryContext& gctx);

  /// Construction path using BarrelEndcapAssembler (wraps
  /// ElementLayerAssembler).
  void constructBarrelEndcap(const Acts::GeometryContext& gctx);

  /// Construction path using ElementLayerAssembler directly. For illustration.
  void constructDirectLayer(const Acts::GeometryContext& gctx);

  /// Construction path using SensorLayerAssembler with groupBy. For
  /// illustration: sensors are collected directly and grouped by walking the
  /// parent chain to find the enclosing layer element.
  void constructDirectLayerGrouped(const Acts::GeometryContext& gctx);

  Config m_cfg;
};

}  // namespace ActsExamples
