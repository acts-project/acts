// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {
namespace Telescope {

/// Global method to build the telescope tracking geometry
///
/// @tparam detector_element_t is the actual type of the detector
/// element, each derivative of a TelescopeDetectorElement can be used
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param positions is the offset w of different layers in the longitudinal
/// direction
/// @param offsets is the offset (u, v) of the layers in the transverse plane
/// @param pSize is the plane size
/// @param thickness is the material thickness of each layer
/// @param binValue indicates which axis the planes normals are parallel to
template <typename detector_element_t>
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename detector_element_t::ContextType& gctx,
    std::vector<std::shared_ptr<detector_element_t>>& detectorStore,
    const std::vector<double>& positions, const std::vector<double>& offsets,
    const std::vector<double>& pSize, double thickness,
    Acts::BinningValue binValue = Acts::BinningValue::binZ) {
  using namespace Acts::UnitLiterals;

  // The rectangle bounds
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>(
      Acts::RectangleBounds(pSize[0] * 0.5, pSize[1] * 0.5));

  // Material of the surfaces
  Acts::Material silicon = Acts::Material::fromMassDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::MaterialSlab matProp(silicon, thickness);
  const auto surfaceMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Construct the rotation
  // This assumes the binValue is binX, binY or binZ. No reset is necessary in
  // case of binZ
  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  if (binValue == Acts::BinningValue::binX) {
    rotation.col(0) = Acts::Vector3(0, 0, -1);
    rotation.col(1) = Acts::Vector3(0, 1, 0);
    rotation.col(2) = Acts::Vector3(1, 0, 0);
  } else if (binValue == Acts::BinningValue::binY) {
    rotation.col(0) = Acts::Vector3(1, 0, 0);
    rotation.col(1) = Acts::Vector3(0, 0, -1);
    rotation.col(2) = Acts::Vector3(0, 1, 0);
  }

  // Construct the surfaces and layers
  size_t nLayers = positions.size();
  std::vector<Acts::LayerPtr> layers(nLayers);
  unsigned int i;
  for (i = 0; i < nLayers; ++i) {
    // The translation without rotation yet
    Acts::Translation3 trans(offsets[0], offsets[1], positions[i]);
    // The transform
    Acts::Transform3 trafo(rotation * trans);
    // Create the detector element
    auto detElement = std::make_shared<TelescopeDetectorElement>(
        Identifier(i), std::make_shared<const Acts::Transform3>(trafo), rBounds,
        1._um, surfaceMaterial);
    // Get the surface
    auto surface = detElement->surface().getSharedPtr();
    // Add the detector element to the event store
    detectorStore.push_back(std::move(detElement));
    // Construct the surface array (one surface contained)
    std::unique_ptr<Acts::SurfaceArray> surArray(
        new Acts::SurfaceArray(surface));
    // Construct the layer
    layers[i] =
        Acts::PlaneLayer::create(trafo, rBounds, std::move(surArray), 1._mm);
    // Associate the layer to the surface
    auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
    mutableSurface->associateLayer(*layers[i]);
  }

  // The volume transform
  Acts::Translation3 transVol(offsets[0], offsets[1],
                              (positions.front() + positions.back()) * 0.5);
  Acts::Transform3 trafoVol(rotation * transVol);

  // The volume bounds is set to be a bit larger than cubic with planes
  auto length = positions.back() - positions.front();
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(
      pSize[0] * 0.5 + 5._mm, pSize[1] * 0.5 + 5._mm, length + 10._mm);

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  Acts::LayerVector layVec;
  for (i = 0; i < nLayers; i++) {
    layVec.push_back(layers[i]);
  }
  // Create the layer array
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      gctx, layVec, positions.front() - 2._mm, positions.back() + 2._mm,
      Acts::BinningType::arbitrary, binValue));

  // Build the tracking volume
  auto trackVolume =
      Acts::TrackingVolume::create(trafoVol, boundsVol, nullptr,
                                   std::move(layArr), nullptr, {}, "Telescope");

  // Build and return tracking geometry
  return std::unique_ptr<Acts::TrackingGeometry>(
      new Acts::TrackingGeometry(trackVolume));
}

}  // end of namespace Telescope
}  // end of namespace ActsExamples
