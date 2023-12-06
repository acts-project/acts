// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/AlignedTelescopeDetector/BuildAlignedTelescopeDetector.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

    /// Gaussian module parameters - 6 Degrees of freedom
double gSigmaX = 0.;  // smear position along local x Axis
double gSigmaY = 0.;  // smear position along local y Axis
double gSigmaZ = 0.;  // smear position along local z Axis
double aSigmaX = 0.;  // rotate around local x Axis
double aSigmaY = 0.;  // rotate around local y Axis
double aSigmaZ = 0.;  // rotate around local z Axis

void applyTransform(Acts::Transform3& trf, std::default_random_engine& rng, double sigmaInPlane, double sigmaOutPlane,
    double sigmaOutRot, double sigmaInRot ) {
    std::normal_distribution<double> gauss(0., 1.);

    // double gSigmaX = 200.;  // smear position along local x Axis
    // double gSigmaY = 200.;  // smear position along local y Axis
    // double gSigmaZ = 0.;  // smear position along local z Axis
    // double aSigmaX = 0.;  // rotate around local x Axis
    // double aSigmaY = 0.;  // rotate around local y Axis
    // double aSigmaZ = 0.54;  // rotate around local z Axis
    gSigmaX = sigmaInPlane;  // smear position along local x Axis
    gSigmaY = sigmaInPlane;  // smear position along local y Axis
    gSigmaZ = sigmaOutPlane;  // smear position along local z Axis
    aSigmaX = sigmaOutRot;  // rotate around local x Axis
    aSigmaY = sigmaOutRot;  // rotate around local y Axis
    aSigmaZ = sigmaInRot;  // rotate around local z Axis
    // the shifts in x, y, z
    double tx = gSigmaX != 0 ? gSigmaX * gauss(rng) : 0.;
    double ty = gSigmaY != 0 ? gSigmaY * gauss(rng) : 0.;
    double tz = gSigmaZ != 0 ? gSigmaZ * gauss(rng) : 0.;

    // Add a translation - if there is any
    if (tx != 0. or ty != 0. or tz != 0.) {
        const auto& tMatrix = trf.matrix();
        auto colX = tMatrix.block<3, 1>(0, 0).transpose();
        auto colY = tMatrix.block<3, 1>(0, 1).transpose();
        auto colZ = tMatrix.block<3, 1>(0, 2).transpose();
        Acts::Vector3 newCenter = tMatrix.block<3, 1>(0, 3).transpose() +
                                  tx * colX + ty * colY + tz * colZ;
        trf.translation() = newCenter;
    }

    // now modify it - rotation around local X
    if (aSigmaX != 0.) {
        trf *= Acts::AngleAxis3(aSigmaX * gauss(rng), Acts::Vector3::UnitX());
    }
    if (aSigmaY != 0.) {
        trf *= Acts::AngleAxis3(aSigmaY * gauss(rng), Acts::Vector3::UnitY());
    }
    if (aSigmaZ != 0.) {
        trf *= Acts::AngleAxis3(aSigmaZ * gauss(rng), Acts::Vector3::UnitZ());
    }
}

std::unique_ptr<const Acts::TrackingGeometry>


ActsExamples::AlignedTelescope::buildDetector(
    const typename ActsExamples::AlignedTelescope::AlignedTelescopeDetectorElement::
        ContextType& gctx,
    std::vector<
        std::shared_ptr<ActsExamples::AlignedTelescope::AlignedTelescopeDetectorElement>>&
        detectorStore,
    const std::vector<double>& positions, const std::array<double, 2>& offsets,
    const std::array<double, 2>& bounds, double thickness,
    ActsExamples::AlignedTelescope::AlignedTelescopeSurfaceType surfaceType,
    Acts::BinningValue binValue, int rnd, double sigmaInPlane, double sigmaOutPlane,
    double sigmaOutRot, double sigmaInRot) {

  using namespace Acts::UnitLiterals;

  // The rectangle bounds for plane surface
  std::default_random_engine rng(rnd);
  std::normal_distribution<double> gauss(0., 1.);
  // const std::array<double, 2> offsets_bis={{offsets[0]*gauss(rng),offsets[1]*gauss(rng)}};

  const auto pBounds =
      std::make_shared<const Acts::RectangleBounds>(bounds[0], bounds[1]);
  // The radial bounds for disc surface
  const auto rBounds =
      std::make_shared<const Acts::RadialBounds>(bounds[0], bounds[1]);

  // Material of the surfaces
  // Acts::Material silicon = Acts::Material::fromMassDensity(
  //     9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3);

  Acts::Material silicon = Acts::Material::fromMassDensity(
      42.20_cm, 90.2_cm, 12.0855, 6, 1.18_g / 1_cm3);

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
  for (unsigned int i = 0; i < nLayers; ++i) {
    // The translation without rotation yet
    // const std::array<double, 2> offsets_bis={{offsets[0]*gauss(rng),offsets[1]*gauss(rng)}};
    // Acts::Translation3 trans(offsets_bis[0], offsets_bis[1], positions[i]);
    // std::default_random_engine rng(315452);
    Acts::Translation3 trans(offsets[0], offsets[1], positions[i]);

    // The transform

    Acts::Transform3 trafo(rotation * trans);

    applyTransform(trafo, rng,sigmaInPlane, sigmaOutPlane, sigmaOutRot, sigmaInRot);
    // Create the detector element
    std::shared_ptr<AlignedTelescopeDetectorElement> detElement = nullptr;
    if (surfaceType == AlignedTelescopeSurfaceType::Plane) {
      detElement = std::make_shared<AlignedTelescopeDetectorElement>(
          std::make_shared<const Acts::Transform3>(trafo), pBounds, 1._um,
          surfaceMaterial);
    } else {
      detElement = std::make_shared<AlignedTelescopeDetectorElement>(
          std::make_shared<const Acts::Transform3>(trafo), rBounds, 1._um,
          surfaceMaterial);
    }
    // Get the surface
    auto surface = detElement->surface().getSharedPtr();
    // Add the detector element to the detector store
    detectorStore.push_back(std::move(detElement));
    // Construct the surface array (one surface contained)
    std::unique_ptr<Acts::SurfaceArray> surArray(
        new Acts::SurfaceArray(surface));
    // Construct the layer
    if (surfaceType == AlignedTelescopeSurfaceType::Plane) {
      layers[i] =
          Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm);
    } else {
      layers[i] =
          Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 1._mm);
    }
    // Associate the layer to the surface
    auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
    mutableSurface->associateLayer(*layers[i]);
  }

  // The volume transform
  Acts::Translation3 transVol(offsets[0], offsets[1],
                              (positions.front() + positions.back()) * 0.5);
  Acts::Transform3 trafoVol(rotation * transVol);

  // The volume bounds is set to be a bit larger than either cubic with planes
  // or cylinder with discs
  auto length = positions.back() - positions.front();
  Acts::VolumeBoundsPtr boundsVol = nullptr;
  if (surfaceType == AlignedTelescopeSurfaceType::Plane) {
    boundsVol = std::make_shared<const Acts::CuboidVolumeBounds>(
        bounds[0] + 5._mm, bounds[1] + 5._mm, length + 10._mm);
  } else {
    boundsVol = std::make_shared<const Acts::CylinderVolumeBounds>(
        std::max(bounds[0] - 5.0_mm, 0.), bounds[1] + 5._mm, length + 10._mm);
  }

  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  Acts::LayerVector layVec;
  for (unsigned int i = 0; i < nLayers; i++) {
    layVec.push_back(layers[i]);
  }
  // Create the layer array
  Acts::GeometryContext genGctx{gctx};
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
      genGctx, layVec, positions.front() - 2._mm, positions.back() + 2._mm,
      Acts::BinningType::arbitrary, binValue));

  // Build the tracking volume
  auto trackVolume =
      Acts::TrackingVolume::create(trafoVol, boundsVol, nullptr,
                                   std::move(layArr), nullptr, {}, "AlignedTelescope");

  // Build and return tracking geometry
  return std::make_unique<Acts::TrackingGeometry>(trackVolume);
}
