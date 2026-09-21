// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ILayerArrayCreator.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/PadBlueprintNode.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Diagnostics.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::Blueprint;
using Acts::BlueprintNode;
using Acts::BlueprintOptions;
using Acts::LayerBlueprintNode;
using Acts::MaterialDesignatorBlueprintNode;
using Acts::PadBlueprintNode;
using Acts::StaticBlueprintNode;

std::unique_ptr<const Acts::TrackingGeometry>
ActsExamples::buildTelescopeDetectorGen3(
    const Acts::GeometryContext& gctx,
    std::vector<std::shared_ptr<const Acts::SurfacePlacementBase>>&
        detectorStore,
    const std::vector<double>& positions,
    const std::vector<double>& stereoAngles,
    const std::array<double, 2>& offsets, const std::array<double, 2>& bounds,
    double thickness, TelescopeSurfaceType surfaceType,
    Acts::AxisDirection rotDirection) {
  using namespace Acts::UnitLiterals;
  std::cout << "Building telescope detector Gen3..." << std::endl;

  auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

  // The rectangle bounds for plane surface
  const auto pBounds =
      std::make_shared<const Acts::RectangleBounds>(bounds[0], bounds[1]);

  // Material of the surfaces
  Acts::Material silicon = Acts::Material::fromMassDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::MaterialSlab matProp(silicon, thickness);
  const auto surfaceMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  //   // temporary
  //   if (rotDirection != Acts::AxisDirection::AxisZ) {
  //     throw std::invalid_argument(
  //         "Only AxisDirection::AxisZ is currently supported, as a possible "
  //         "rotation");
  //   }

  // Construct the rotation
  // This assumes the direction is AxisX, AxisY or AxisZ. No reset is necessary
  // in case of AxisZ
  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  if (rotDirection == Acts::AxisDirection::AxisX) {
    rotation.col(0) = Acts::Vector3(0, 0, -1);
    rotation.col(1) = Acts::Vector3(0, 1, 0);
    rotation.col(2) = Acts::Vector3(1, 0, 0);
  } else if (rotDirection == Acts::AxisDirection::AxisY) {
    rotation.col(0) = Acts::Vector3(1, 0, 0);
    rotation.col(1) = Acts::Vector3(0, 0, -1);
    rotation.col(2) = Acts::Vector3(0, 1, 0);
  }

  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisX] = {200_mm, 200_mm};
  cfg.envelope[AxisDirection::AxisY] = {200_mm, 200_mm};
  cfg.envelope[AxisDirection::AxisZ] = {200_mm, 200_mm};
  // cfg.envelope[AxisDirection::AxisR] = {1_mm, 2_mm};
  Blueprint root{cfg};

  auto& cubcontainer = root.addCuboidContainer("CuboidContainer", rotDirection);

  std::filesystem::path debugOutputDir = "telescope_debug_obj";
  std::filesystem::create_directories(debugOutputDir);

  std::size_t nLayers = positions.size();
  for (unsigned int i = 0; i < nLayers; i++) {
    // The translation without rotation yet
    Acts::Vector3 transVec;
    if (rotDirection == Acts::AxisDirection::AxisX) {
      transVec = Acts::Vector3(positions[i], offsets[0], offsets[1]);
    } else if (rotDirection == Acts::AxisDirection::AxisY) {
      transVec = Acts::Vector3(offsets[0], positions[i], offsets[1]);
    } else {  // AxisZ
      transVec = Acts::Vector3(offsets[0], offsets[1], positions[i]);
    }
    Acts::Translation3 trans(transVec);

    // The entire transformation (the coordinate system, whose center is defined
    // by trans, will be rotated as well)
    Acts::Transform3 trafo(trans * rotation);
    // Acts::Transform3 trafo(rotation * trans);  // this is the wrong order of
    // matrix multiplication

    // rotate around local z axis by stereo angle
    auto stereo = stereoAngles[i];
    trafo *= Acts::AngleAxis3(stereo, Acts::Vector3::UnitZ());

    // std::cout << trafo.matrix() << std::endl;

    // Create the detector element
    std::shared_ptr<TelescopeDetectorElement> detElement = nullptr;

    const auto id =
        static_cast<TelescopeDetectorElement::Identifier>(detectorStore.size());

    if (surfaceType == TelescopeSurfaceType::Plane) {
      detElement = std::make_shared<TelescopeDetectorElement>(
          id, std::make_shared<Acts::Transform3>(trafo), pBounds, 1._um,
          surfaceMaterial);
    } else {
      throw std::invalid_argument(
          "Only TelescopeSurfaceType::Plane is currently supported, not "
          "TelescopeSurfaceType::Disc");
      //   detElement = std::make_shared<TelescopeDetectorElement>(
      //       id, std::make_shared<Acts::Transform3>(trafo), rBounds, 1._um,
      //       surfaceMaterial);
    }
    detectorStore.push_back(detElement);

    auto layerBounds = std::make_shared<CuboidVolumeBounds>(
        bounds[0] + 5_mm, bounds[1] + 5_mm, thickness / 2 + 1_mm);

    if (rotDirection == Acts::AxisDirection::AxisX) {
      layerBounds = std::make_shared<CuboidVolumeBounds>(
          thickness / 2 + 1_mm, bounds[0] + 5_mm, bounds[1] + 5_mm);
    } else if (rotDirection == Acts::AxisDirection::AxisY) {
      layerBounds = std::make_shared<CuboidVolumeBounds>(
          bounds[0] + 5_mm, thickness / 2 + 1_mm, bounds[1] + 5_mm);
    } else if (rotDirection == Acts::AxisDirection::AxisZ) {
      layerBounds = std::make_shared<CuboidVolumeBounds>(
          bounds[0] + 5_mm, bounds[1] + 5_mm, thickness / 2 + 1_mm);
    }

    auto layerVol = std::make_unique<TrackingVolume>(
        Acts::Transform3{trans}, layerBounds, "parent" + std::to_string(i));

    // Get the surface
    auto surface = detElement->surface().getSharedPtr();
    layerVol->addSurface(surface);

    // // --- Debug-Export ---
    // {
    //   Acts::ObjVisualization3D vis;
    //   Acts::GeometryView3D::drawVolume(vis, *layerVol, gctx);
    //   vis.write((debugOutputDir / ("layer_" + std::to_string(i) +
    //   ".obj")).string()); Acts::GeometryView3D::drawSurface(vis, *surface,
    //   gctx); vis.write((debugOutputDir / ("surface_" + std::to_string(i) +
    //   ".obj")).string());
    // }
    // // --- End Debug-Export ---

    auto layerNode = std::make_shared<StaticBlueprintNode>(std::move(layerVol));

    cubcontainer.addChild(std::move(layerNode));
  }

  //   std::ofstream os{"telescope.dot"};
  //   root.graphviz(os);
  auto trackingGeometry = root.construct({}, gctx, *logger);

  return trackingGeometry;
}
