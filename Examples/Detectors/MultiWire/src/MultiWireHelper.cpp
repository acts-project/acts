// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MultiWire/MultiWireHelper.hpp"

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/VolumeStructureBuilder.hpp"
#include "Acts/Detector/interface/IExternalStructureBuilder.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "ActsExamples/Geant4/GdmlDetectorConstruction.hpp"
#include "ActsExamples/Geant4Detector/Geant4Detector.hpp"

#include <iostream>
#include <vector>

std::vector<std::shared_ptr<Acts::Surface>>
ActsExamples::MultiWireHelper::getStrawSurfaces(
    std::vector<std::string> sensitiveNames,
    std::vector<std::string> passiveNames) {
  const Acts::GeometryContext gctx;

  std::string gdmlPath =
      "../../acts/Examples/Detectors/MuonSpectrometerMockupDetector/"
      "MuonChamber.gdml";

  // Geant4Detector Config creator with the g4world from the gdml file
  ActsExamples::GdmlDetectorConstruction geo_gdml(gdmlPath);
  auto g4WorldConfig = ActsExamples::Geant4::Geant4Detector::Config();
  g4WorldConfig.name = "MultiLayerChamber";
  g4WorldConfig.g4World = geo_gdml.Construct();

  // Get the sensitive and passive surfaces and pass to the g4World Config
  auto g4Sensitive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          sensitiveNames);
  auto g4Passive =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          passiveNames);

  auto g4SurfaceOptions = Acts::Geant4DetectorSurfaceFactory::Options();
  g4SurfaceOptions.sensitiveSurfaceSelector = g4Sensitive;
  g4SurfaceOptions.passiveSurfaceSelector = g4Passive;
  g4WorldConfig.g4SurfaceOptions = g4SurfaceOptions;

  auto g4detector = ActsExamples::Geant4::Geant4Detector();

  auto [detector, surfaces, detectorElements] =
      g4detector.constructDetector(g4WorldConfig, Acts::getDummyLogger());

  std::vector<std::shared_ptr<Acts::Surface>> strawSurfaces;

  strawSurfaces.reserve(detectorElements.size());

  // Convert the physical volumes of the detector elements to straw surfaces
  std::cout << detectorElements.size() << std::endl;
  std::cin.ignore();
  for (auto& detectorElement : detectorElements) {
    // auto context = Acts::GeometryContext();
    auto g4conv = Acts::Geant4PhysicalVolumeConverter();

    g4conv.forcedType = Acts::Surface::SurfaceType::Straw;
    auto g4ConvSurf = g4conv.Geant4PhysicalVolumeConverter::surface(
        detectorElement->g4PhysicalVolume(), detectorElement->transform(gctx));

    strawSurfaces.push_back(g4ConvSurf);
  }

  return strawSurfaces;
}

std::vector<Acts::Experimental::LayerStructureBuilder::Binning>
ActsExamples::MultiWireHelper::layerBinning(
    std::vector<std::shared_ptr<Acts::Surface>> surfaces,
    std::array<std::pair<float, float>, 3> multiWireBounds) {
  const Acts::GeometryContext context;

  // Find the number of surfaces along y and z
  // A surface as reference
  auto refSurf = surfaces.front();
  auto tolerance = surfaces.front()->bounds().values()[0];
  std::size_t nSurfacesZ = 0, nSurfacesY = 0;

  for (auto& surf : surfaces) {
    if (abs(surf->center(context).y() - refSurf->center(context).y()) <=
        tolerance) {
      nSurfacesZ += 1;
    };
    if (abs(surf->center(context).z() - refSurf->center(context).z()) <=
        tolerance) {
      nSurfacesY += 1;
    };
  }

  using LayerBinning = Acts::Experimental::LayerStructureBuilder::Binning;

  std::vector<LayerBinning> lbinnings = {
      LayerBinning{Acts::BinningData(Acts::closed, Acts::binZ, nSurfacesZ,
                                     multiWireBounds[2].first,
                                     multiWireBounds[2].second),
                   1u},
      LayerBinning{Acts::BinningData(Acts::closed, Acts::binY, nSurfacesY,
                                     multiWireBounds[1].first,
                                     multiWireBounds[1].second),
                   0u}};

  return lbinnings;
}

std::shared_ptr<Acts::Experimental::LayerStructureBuilder>
ActsExamples::MultiWireHelper::internalLayerBuilder(
    Acts::Experimental::LayerStructureBuilder::Config lConfig) {
  auto iLayerBuilder =
      std::make_shared<Acts::Experimental::LayerStructureBuilder>(
          lConfig, Acts::getDefaultLogger("InternalLayerBuilder",
                                          Acts::Logging::VERBOSE));
  return iLayerBuilder;
}

std::shared_ptr<Acts::Experimental::VolumeStructureBuilder>
ActsExamples::MultiWireHelper::externalVolumeBuilder(
    Acts::Experimental::VolumeStructureBuilder::Config vsConfig) {
  auto vBuilder = std::make_shared<Acts::Experimental::VolumeStructureBuilder>(
      vsConfig,
      Acts::getDefaultLogger("ExternalVolumeBuilder", Acts::Logging::VERBOSE));
  return vBuilder;
}

std::array<std::pair<float, float>, 3>
ActsExamples::MultiWireHelper::getMultiWireBounds(
    std::vector<std::shared_ptr<Acts::Surface>> surfaces) {
  Acts::GeometryContext context;

  std::array<std::pair<float, float>, 3> min_max;
  std::fill(min_max.begin(), min_max.end(),
            std::make_pair<float, float>(std::numeric_limits<float>::max(),
                                         -std::numeric_limits<float>::max()));
  for (auto& surf : surfaces) {
    min_max[0].first =
        std::min(min_max[0].first, (float)surf->center(context).x());
    min_max[0].second =
        std::max(min_max[0].second, (float)surf->center(context).x());

    min_max[1].first =
        std::min(min_max[1].first, (float)surf->center(context).y());
    min_max[1].second =
        std::max(min_max[1].second, (float)surf->center(context).y());

    min_max[2].first =
        std::min(min_max[2].first, (float)surf->center(context).z());
    min_max[2].second =
        std::max(min_max[2].second, (float)surf->center(context).z());
  }

  return min_max;
}
