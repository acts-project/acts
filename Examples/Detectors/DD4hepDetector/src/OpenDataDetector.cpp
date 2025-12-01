// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
#include <Acts/Geometry/Blueprint.hpp>
#include <Acts/Geometry/BlueprintOptions.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp>

#include <algorithm>
#include <iterator>

#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>

namespace ActsExamples {

OpenDataDetector::OpenDataDetector(const Config& cfg,
                                   const Acts::GeometryContext& gctx)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  ACTS_INFO("OpenDataDetector construct");
  construct(gctx);
}

auto OpenDataDetector::config() const -> const Config& {
  return m_cfg;
}

std::shared_ptr<ActsPlugins::DD4hepDetectorElement>
OpenDataDetector::defaultDetectorElementFactory(
    const dd4hep::DetElement& element, const std::string& axes, double scale) {
  return std::make_shared<ActsPlugins::DD4hepDetectorElement>(element, axes,
                                                              scale);
}

void OpenDataDetector::construct(const Acts::GeometryContext& gctx) {
  using namespace Acts::Experimental;
  using namespace Acts;
  using namespace Acts::UnitLiterals;
  using enum AxisDirection;

  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {
          .dd4hepDetector = &dd4hepDetector(),
          .lengthScale = Acts::UnitConstants::cm,
      },
      logger().cloneWithSuffix("BlpBld")};

  // func(world);

  // ACTS_VERBOSE("Found detelement: " << pixelBrl0.name());

  // std::vector sensitives = resolveSensitives(pixelBrl0);

  // std::vector<Acts::Surface*> surfaces;
  // surfaces.reserve(sensitives.size());
  // std::vector<std::shared_ptr<ActsPlugins::DD4hepDetectorElement>>
  //     detectorElements;
  // BARREL: XYZ
  // ENDCAP: XZY
  // for (const auto& sensitive : sensitives) {
  //   builder.createDetectorElement(sensitive, "XYZ");
  // }

  Blueprint::Config cfg;
  cfg.envelope[AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisR] = {0_mm, 20_mm};
  Blueprint root{cfg};

  dd4hep::DetElement pixelBarrelElem =
      builder.findDetElementByName("PixelBarrel").value();

  auto pixelBarrel = builder.addLayers(pixelBarrelElem, "XYZ", AxisR,
                                       std::regex{"PixelLayer\\d"});

  root.addChild(pixelBarrel);

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
