// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/OpenDataDetector.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include <Acts/Geometry/Blueprint.hpp>
#include <Acts/Geometry/BlueprintOptions.hpp>
#include <Acts/Geometry/BoundarySurfaceFace.hpp>
#include <Acts/Geometry/ContainerBlueprintNode.hpp>
#include <Acts/Geometry/CylinderVolumeBounds.hpp>
#include <Acts/Geometry/Extent.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Navigation/CylinderNavigationPolicy.hpp>
#include <Acts/Navigation/SurfaceArrayNavigationPolicy.hpp>
#include <Acts/Navigation/TryAllNavigationPolicy.hpp>
#include <Acts/Surfaces/SurfaceArray.hpp>
#include <Acts/Utilities/AxisDefinitions.hpp>
#include <ActsPlugins/DD4hep/DD4hepConversionHelpers.hpp>

#include <algorithm>
#include <iterator>
#include <utility>

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

namespace {

int parseLayerNumber(const dd4hep::DetElement& elem,
                     const std::regex& pattern) {
  std::cmatch match;

  if (!std::regex_match(elem.name(), match, pattern)) {
    throw std::runtime_error(std::format(
        "Layer name {} does not match expected pattern", elem.name()));
  }

  if (match.size() != 2) {
    throw std::runtime_error(std::format(
        "Layer name {} matched pattern but did not capture layer number",
        elem.name()));
  }

  int n = std::stoi(match[1]);
  return n;
}

}  // namespace

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

  // BARREL: XYZ
  // ENDCAP: XZY

  Blueprint::Config cfg;
  cfg.envelope[AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisR] = {0_mm, 20_mm};
  Blueprint root{cfg};

  auto& outer = root.addCylinderContainer("OpenDataDetector", AxisR);

  outer.addChild(builder.makeBeampipe());

  using AttachmentStrategy = Acts::VolumeAttachmentStrategy;
  using SrfArrayNavPol = Acts::SurfaceArrayNavigationPolicy;

  auto constant = [this](const std::string& name) -> int {
    return dd4hepDetector().constant<int>(name);
  };

  auto makeCustomizer = [&](const std::string& det,
                            const std::regex& layerPattern)
      -> ActsPlugins::DD4hep::LayerHelper::Customizer {
    return [&, det](const auto& elem,
                    Acts::Experimental::LayerBlueprintNode& layer) {
      int n = parseLayerNumber(elem, layerPattern);

      std::string name = elem.name();
      using enum SrfArrayNavPol::LayerType;
      SrfArrayNavPol::Config cfg;

      if (name.find("Endcap") != std::string::npos) {
        cfg.bins = {constant(std::format("{}_e_sf_b_r", det)),
                    constant(std::format("{}_e_sf_b_phi", det))};
        cfg.layerType = Disc;
      } else {
        cfg.bins = {constant(std::format("{}_b{}_sf_b_phi", det, n)),
                    constant(std::format("{}_b_sf_b_z", det))};
        cfg.layerType = Cylinder;
      }

      layer.setNavigationPolicyFactory(NavigationPolicyFactory{}
                                           .add<CylinderNavigationPolicy>()
                                           .add<SrfArrayNavPol>(cfg)
                                           .asUniquePtr());

      // @TODO: This needs to set the material as well
    };
  };

  auto pixelAssembly = builder.findDetElementByName("Pixels").value();
  std::regex pixelLayerPattern{"(?:PixelLayer|PixelEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly(pixelAssembly)
      .setBarrelAxes("XYZ")
      .setEndcapAxes("XZY")
      .setLayerPattern(pixelLayerPattern)
      .setCustomizer(makeCustomizer("pix", pixelLayerPattern))
      .addTo(outer);

  auto sstripAssembly = builder.findDetElementByName("ShortStrips").value();
  std::regex sstripLayerPattern{
      "(?:ShortStripLayer|ShortStripEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly(sstripAssembly)
      .setBarrelAxes("XYZ")
      .setEndcapAxes("XZY")
      .setLayerPattern(sstripLayerPattern)
      .setCustomizer(makeCustomizer("ss", sstripLayerPattern))
      .addTo(outer);

  auto lstripAssembly = builder.findDetElementByName("LongStrips").value();
  std::regex lstripLayerPattern{"(?:LongStripLayer|LongStripEndcap[NP])(\\d)"};
  builder.barrelEndcapAssemblyHelper()
      .setAssembly(lstripAssembly)
      .setBarrelAxes("XYZ")
      .setEndcapAxes("XZY")
      .setLayerPattern(lstripLayerPattern)
      .setCustomizer(makeCustomizer("ls", lstripLayerPattern))
      .addTo(outer);

  BlueprintOptions options;

  m_trackingGeometry = root.construct(options, gctx, logger());
}

}  // namespace ActsExamples
