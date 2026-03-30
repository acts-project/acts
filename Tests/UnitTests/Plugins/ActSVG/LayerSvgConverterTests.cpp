// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/ActSVG/LayerSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <fstream>
#include <memory>
#include <numbers>
#include <vector>

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

auto tgContext = GeometryContext::dangerouslyDefaultConstruct();

std::shared_ptr<const LayerCreator> lCreator(nullptr);

void setupTools() {
  if (lCreator == nullptr) {
    LayerCreator::Config lCreatorCfg;
    lCreatorCfg.surfaceArrayCreator =
        std::make_shared<const SurfaceArrayCreator>();
    lCreator = std::make_shared<const LayerCreator>(lCreatorCfg);
  }
}

std::shared_ptr<Layer> generateDiscLayer(double rInner, double rOuter,
                                         unsigned int quarterSegments,
                                         unsigned int nRings,
                                         bool useTrapezoids = false) {
  // Some preparations
  setupTools();
  unsigned int fullSegments = 4 * quarterSegments;
  std::vector<std::shared_ptr<const Surface>> moduleSurfaces;
  double phiStep = 2 * std::numbers::pi / fullSegments;
  double rStep = (rOuter - rInner) / nRings;
  // Reserve & fill
  moduleSurfaces.reserve(fullSegments * nRings);
  // Radial disc
  if (!useTrapezoids) {
    for (unsigned int ir = 0; ir < nRings; ++ir) {
      std::shared_ptr<const RadialBounds> rBounds = nullptr;
      rBounds = std::make_shared<RadialBounds>(
          rInner + ir * rStep - 0.025 * rInner,
          rInner + (ir + 1u) * rStep + 0.025 * rInner, 0.55 * phiStep, 0.);
      for (unsigned int is = 0; is < fullSegments; ++is) {
        // Place the module
        auto placement = Transform3::Identity();
        if ((is % 2) != 0u) {
          placement.pretranslate(Vector3{0., 0., 2.});
        }
        placement.rotate(Eigen::AngleAxisd(is * phiStep, Vector3(0, 0, 1)));
        auto dModule = Surface::makeShared<DiscSurface>(placement, rBounds);
        moduleSurfaces.push_back(dModule);
      }
    }
  } else {
    for (unsigned int ir = 0; ir < nRings; ++ir) {
      // Trapezoid parameters
      double radius = rInner + (ir + 0.5) * rStep;
      double yHalf = rStep * 0.5125;

      double xHalfMin =
          1.15 * (rInner + ir * rStep) * std::numbers::pi / fullSegments;
      double xHalfMax =
          1.15 * (rInner + (ir + 1) * rStep) * std::numbers::pi / fullSegments;

      std::shared_ptr<const TrapezoidBounds> tBounds =
          std::make_shared<const TrapezoidBounds>(xHalfMin, xHalfMax, yHalf);
      for (unsigned int is = 0; is < fullSegments; ++is) {
        // Setting the phi
        double cphi = -std::numbers::pi + is * phiStep;
        Vector3 center(radius * std::cos(cphi), radius * std::sin(cphi),
                       (is % 2) * 2 + (ir % 2) * 5);
        // Local axis system
        Vector3 localY(std::cos(cphi), std::sin(cphi), 0.);
        Vector3 localZ(0., 0., 1.);
        Vector3 localX = localY.cross(localZ);
        RotationMatrix3 rotation;
        rotation.col(0) = localX;
        rotation.col(1) = localY;
        rotation.col(2) = localZ;
        Transform3 placement(Translation3(center) * rotation);
        // Create the module surface
        auto dModule = Surface::makeShared<PlaneSurface>(placement, tBounds);
        moduleSurfaces.push_back(dModule);
      }
    }
  }
  // Let's create the disc layer
  return lCreator->discLayer(tgContext, moduleSurfaces, nRings, fullSegments);
}

BOOST_AUTO_TEST_SUITE(ActSvgSuite)

BOOST_AUTO_TEST_CASE(DiscLayerRadialSvg) {
  // Planar style
  Svg::Style discLayerStyle;
  discLayerStyle.fillColor = {51, 153, 255};
  discLayerStyle.fillOpacity = 0.75;
  discLayerStyle.highlightColor = {255, 153, 51};
  discLayerStyle.highlights = {"mouseover", "mouseout"};
  discLayerStyle.strokeColor = {25, 25, 25};
  discLayerStyle.strokeWidth = 0.5;
  discLayerStyle.quarterSegments = 72u;

  GeometryIdentifier geoID{0};

  // Get the layer
  auto discLayer = generateDiscLayer(100, 250, 32u, 4u);

  Svg::LayerConverter::Options lOptions;
  lOptions.name = "disc_layer_sectors";
  lOptions.surfaceStyles =
      GeometryHierarchyMap<Svg::Style>({{geoID, discLayerStyle}});

  // Get the layer sheets
  auto discLayerSheets =
      Svg::LayerConverter::convert(tgContext, *discLayer, lOptions);

  for (const auto& s : discLayerSheets) {
    Svg::toFile({s}, s._id + ".svg");
  }
}

BOOST_AUTO_TEST_CASE(DiscLayerTrapezoidSvg) {
  // Planar style
  Svg::Style discLayerStyle;
  discLayerStyle.fillColor = {51, 153, 255};
  discLayerStyle.fillOpacity = 0.75;
  discLayerStyle.highlightColor = {255, 153, 51};
  discLayerStyle.highlights = {"mouseover", "mouseout"};
  discLayerStyle.strokeColor = {25, 25, 25};
  discLayerStyle.strokeWidth = 0.5;
  discLayerStyle.quarterSegments = 72u;

  GeometryIdentifier geoID{0};

  // Get the layer
  auto discLayer = generateDiscLayer(100, 250, 32u, 4u, true);

  Svg::LayerConverter::Options lOptions;
  lOptions.name = "disc_layer_trapezoid";
  lOptions.surfaceStyles =
      GeometryHierarchyMap<Svg::Style>({{geoID, discLayerStyle}});

  // Get the layer sheets
  auto discLayerSheets =
      Svg::LayerConverter::convert(tgContext, *discLayer, lOptions);

  for (const auto& s : discLayerSheets) {
    Svg::toFile({s}, s._id + ".svg");
  }
}

BOOST_AUTO_TEST_CASE(CylinderLayerSvg) {
  // Planar style
  Svg::Style cylinderLayerStyle;
  cylinderLayerStyle.fillColor = {51, 153, 255};
  cylinderLayerStyle.fillOpacity = 0.75;
  cylinderLayerStyle.highlightColor = {255, 153, 51};
  cylinderLayerStyle.highlights = {"mouseover", "mouseout"};
  cylinderLayerStyle.strokeColor = {25, 25, 25};
  cylinderLayerStyle.strokeWidth = 0.5;
  cylinderLayerStyle.quarterSegments = 72u;

  GeometryIdentifier geoID{0};

  CylindricalTrackingGeometry cGeometry(tgContext);
  auto tGeometry = cGeometry();
  auto pixelVolume =
      tGeometry->lowestTrackingVolume(tgContext, Vector3(50., 0., 0.));
  if (pixelVolume != nullptr && pixelVolume->confinedLayers() != nullptr) {
    auto layers = pixelVolume->confinedLayers()->arrayObjects();
    std::size_t il = 0;
    for (const auto& layer : layers) {
      if (layer->surfaceArray() == nullptr) {
        continue;
      }

      Svg::LayerConverter::Options lOptions;
      lOptions.name = "cylinder_layer_" + std::to_string(il++);
      lOptions.surfaceStyles =
          GeometryHierarchyMap<Svg::Style>({{geoID, cylinderLayerStyle}});

      // Get the layer sheets
      auto layerSheets =
          Svg::LayerConverter::convert(tgContext, *layer, lOptions);
      for (const auto& s : layerSheets) {
        Svg::toFile({s}, s._id + ".svg");
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(PlaeyLayerSvg) {}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
