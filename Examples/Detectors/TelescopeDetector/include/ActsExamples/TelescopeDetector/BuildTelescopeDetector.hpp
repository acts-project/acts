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
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
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
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
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
/// @param matDecorator is an optional decorator for the material
template <typename detector_element_t>
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename detector_element_t::ContextType& gctx,
    std::vector<std::shared_ptr<detector_element_t>>& detectorStore,
    const std::vector<double>& layerRelDists,
    const Acts::Vector3& firstLayerCenter, const std::vector<double>& boundary,
    double thickness, Acts::BinningValue binValue = Acts::BinningValue::binX) {
  using namespace Acts::UnitLiterals;

  // The rectangle bounds
  // @todo make the bound size configurable
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>(
      Acts::RectangleBounds(boundary[0], boundary[1]));

  // Material of the surfaces
  Acts::Material silicon = Acts::Material::fromMolarDensity(
      9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);
  Acts::MaterialSlab matProp(silicon, thickness);
  const auto surfaceMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Construct the rotation
  // This assumes the binValue is binX, binY or binZ
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

  // Set translation vectors
  size_t nLayers = layerRelDists.size() + 1;
  std::vector<Acts::Vector3> translations;
  translations.reserve(nLayers);
  translations.push_back(firstLayerCenter);
  for (const auto dist : layerRelDists) {
    Acts::Vector3 shift(0, 0, 0);
    // This assumes the binValue is binX, binY or binZ
    shift[binValue] = dist;
    translations.push_back(firstLayerCenter + shift);
  }

  // Construct surfaces
  std::vector<std::shared_ptr<const Acts::Surface>> surfaces;
  surfaces.reserve(nLayers);
  unsigned int i;
  for (i = 0; i < nLayers; i++) {
    Acts::Transform3 trafo(Acts::Transform3::Identity() * rotation);
    trafo.translation() = translations[i];
    // Create the detector element
    auto detElement = std::make_shared<TelescopeDetectorElement>(
        Identifier(i), std::make_shared<const Acts::Transform3>(trafo), rBounds,
        1_um, surfaceMaterial);
    // And remember the surface
    surfaces.push_back(detElement->surface().getSharedPtr());
    // Add it to the event store
    detectorStore.push_back(std::move(detElement));
  }

  // Construct layers
  std::vector<Acts::LayerPtr> layers;
  layers.reserve(nLayers);
  for (i = 0; i < nLayers; i++) {
    Acts::Transform3 trafo(Acts::Transform3::Identity() * rotation);
    trafo.translation() = translations[i];

    std::unique_ptr<Acts::SurfaceArray> surArray(
        new Acts::SurfaceArray(surfaces[i]));

    layers.push_back(
        Acts::PlaneLayer::create(trafo, rBounds, std::move(surArray), 1._mm));

    auto mutableSurface = const_cast<Acts::Surface*>(surfaces[i].get());
    mutableSurface->associateLayer(*layers[i]);
  }

  // Build volume for the layers
  Acts::Transform3 trafoVol(Acts::Transform3::Identity() * rotation);
  trafoVol.translation() = (translations.front() + translations.back()) * 0.5;

  // The volume bounds is a bit larger than the cubic with layers
  auto boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(
      boundary[0] + 5._mm, boundary[1] + 5._mm, layerRelDists.back() + 10._mm);

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));

  Acts::LayerVector layVec;
  for (i = 0; i < nLayers; i++) {
    layVec.push_back(layers[i]);
  }
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      gctx, layVec, (translations.front())[binValue] - 2._mm,
      (translations.back())[binValue] + 2._mm, Acts::BinningType::arbitrary,
      binValue));

  auto trackVolume =
      Acts::TrackingVolume::create(trafoVol, boundsVol, nullptr,
                                   std::move(layArr), nullptr, {}, "Telescope");

  // Build and return tracking geometry
  return std::unique_ptr<Acts::TrackingGeometry>(
      new Acts::TrackingGeometry(trackVolume));
}

}  // end of namespace Telescope
}  // end of namespace ActsExamples
