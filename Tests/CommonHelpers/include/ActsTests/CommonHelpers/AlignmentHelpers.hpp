// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/MeasurementsCreator.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <random>
#include <string>

namespace ActsTests::AlignmentUtils {
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitterType =
    Acts::KalmanFitter<ConstantFieldPropagator, Acts::VectorMultiTrajectory>;

/// helper struct packaging commonly used data members for alignment
/// unit tests
struct aliTestUtils {
  KalmanUpdater kfUpdater;
  KalmanSmoother kfSmoother;
  const Acts::GeometryContext geoCtx =
      Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::MagneticFieldContext magCtx;
  const Acts::CalibrationContext calCtx;
  std::normal_distribution<double> normalDist{0., 1.};
  std::default_random_engine rng{42};
  // detector resolutions
  const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                          {30_um, 50_um}};
  const MeasurementResolutionMap resolutions = {
      {Acts::GeometryIdentifier(), resPixel},
  };
};

// Create a test context

Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> getExtensions(
    aliTestUtils& utils) {
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  extensions.calibrator.connect<&Acts::detail::Test::testSourceLinkCalibrator<
      Acts::VectorMultiTrajectory>>();
  extensions.updater
      .connect<&KalmanUpdater::operator()<Acts::VectorMultiTrajectory>>(
          &utils.kfUpdater);
  extensions.smoother
      .connect<&KalmanSmoother::operator()<Acts::VectorMultiTrajectory>>(
          &utils.kfSmoother);
  return extensions;
}

///
/// @brief Construct a telescope-like detector
///
struct TelescopeDetector {
  /// Default constructor for the Cubic tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  explicit TelescopeDetector(
      std::reference_wrapper<const Acts::GeometryContext> gctx)
      : geoContext(gctx) {
    // Construct the rotation
    rotation.col(0) = Acts::Vector3(0, 0, -1);
    rotation.col(1) = Acts::Vector3(0, 1, 0);
    rotation.col(2) = Acts::Vector3(1, 0, 0);

    // Boundaries of the surfaces
    rBounds = std::make_shared<const Acts::RectangleBounds>(0.1_m, 0.1_m);

    // Material of the surfaces
    Acts::MaterialSlab matProp(makeSilicon(), 80_um);

    surfaceMaterial =
        std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
  }

  ///
  /// Call operator to build the standard cubic tracking geometry
  ///
  std::shared_ptr<const Acts::TrackingGeometry> operator()() {
    using namespace Acts::UnitLiterals;

    unsigned int nLayers = 6;
    std::vector<double> positions = {-500_mm, -300_mm, -100_mm,
                                     100_mm,  300_mm,  500_mm};
    auto length = positions.back() - positions.front();

    std::vector<Acts::LayerPtr> layers(nLayers);
    for (unsigned int i = 0; i < nLayers; ++i) {
      // The transform
      Acts::Translation3 trans(0., 0., positions[i]);
      Acts::Transform3 trafo(rotation * trans);
      auto detElement = std::make_shared<DetectorElementStub>(
          trafo, rBounds, 1._um, surfaceMaterial);
      // The surface is not right!!!
      auto surface = detElement->surface().getSharedPtr();
      // Add it to the event store
      detectorStore.push_back(std::move(detElement));
      std::unique_ptr<Acts::SurfaceArray> surArray(
          new Acts::SurfaceArray(surface));
      // The layer thickness should not be too large
      layers[i] = Acts::PlaneLayer::create(
          trafo, rBounds, std::move(surArray),
          1._mm);  // Associate the layer to the surface
      auto mutableSurface = const_cast<Acts::Surface*>(surface.get());
      mutableSurface->associateLayer(*layers[i]);
    }

    // The volume transform
    Acts::Translation3 transVol(0, 0, 0);
    Acts::Transform3 trafoVol(rotation * transVol);
    auto boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
        rBounds->halfLengthX() + 10._mm, rBounds->halfLengthY() + 10._mm,
        length + 10._mm);

    Acts::LayerArrayCreator::Config lacConfig;
    Acts::LayerArrayCreator layArrCreator(
        lacConfig,
        Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
    Acts::LayerVector layVec;
    for (unsigned int i = 0; i < nLayers; i++) {
      layVec.push_back(layers[i]);
    }

    // Create the layer array
    std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(
        geoContext, layVec, positions.front() - 2._mm, positions.back() + 2._mm,
        Acts::BinningType::arbitrary, Acts::AxisDirection::AxisX));

    // Build the tracking volume
    auto trackVolume = std::make_shared<Acts::TrackingVolume>(
        trafoVol, boundsVol, nullptr, std::move(layArr), nullptr,
        Acts::MutableTrackingVolumeVector{}, "Telescope");

    return std::make_shared<const Acts::TrackingGeometry>(trackVolume);
  }

  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  std::shared_ptr<const Acts::RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const Acts::ISurfaceMaterial> surfaceMaterial = nullptr;

  std::vector<std::shared_ptr<DetectorElementStub>> detectorStore;

  std::reference_wrapper<const Acts::GeometryContext> geoContext;
};

// Construct a straight-line propagator.
StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{std::move(geo)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::StraightLineStepper stepper;
  return StraightPropagator(stepper, std::move(navigator));
}

// Construct a propagator using a constant magnetic field along z.
ConstantFieldPropagator makeConstantFieldPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo, double bz,
    std::unique_ptr<const Acts::Logger> logger) {
  Acts::Navigator::Config cfg{std::move(geo)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg, logger->cloneWithSuffix("Nav"));
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  ConstantFieldStepper stepper(std::move(field));
  return ConstantFieldPropagator(std::move(stepper), std::move(navigator),
                                 logger->cloneWithSuffix("Prop"));
}

// Construct initial track parameters.
Acts::BoundTrackParameters makeParameters(aliTestUtils& utils) {
  using namespace Acts::UnitLiterals;
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 0.5_degree;
  stddev[Acts::eBoundTheta] = 0.5_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();

  auto loc0 = 0. + stddev[Acts::eBoundLoc0] * utils.normalDist(utils.rng);
  auto loc1 = 0. + stddev[Acts::eBoundLoc1] * utils.normalDist(utils.rng);
  auto t = 42_ns + stddev[Acts::eBoundTime] * utils.normalDist(utils.rng);
  auto phi = 0_degree + stddev[Acts::eBoundPhi] * utils.normalDist(utils.rng);
  auto theta =
      90_degree + stddev[Acts::eBoundTheta] * utils.normalDist(utils.rng);
  auto qOverP =
      1_e / 1_GeV + stddev[Acts::eBoundQOverP] * utils.normalDist(utils.rng);

  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-1_m, loc0, loc1, t);

  return Acts::BoundTrackParameters::createCurvilinear(
      mPos4, phi, theta, qOverP, cov, Acts::ParticleHypothesis::pion());
}

struct KalmanFitterInputTrajectory {
  // The source links
  std::vector<Acts::detail::Test::TestSourceLink> sourceLinks;
  // The start parameters
  std::optional<Acts::BoundTrackParameters> startParameters;
};

///
/// Function to create trajectories for kalman fitter
///
std::vector<KalmanFitterInputTrajectory> createTrajectories(
    std::shared_ptr<const Acts::TrackingGeometry> geo,
    std::size_t nTrajectories, aliTestUtils& util) {
  // simulation propagator
  const auto simPropagator = makeStraightPropagator(std::move(geo));

  std::vector<KalmanFitterInputTrajectory> trajectories;
  trajectories.reserve(nTrajectories);

  for (unsigned int iTrack = 0; iTrack < nTrajectories; iTrack++) {
    auto start = makeParameters(util);
    // Launch and collect - the measurements
    auto measurements =
        createMeasurements(simPropagator, util.geoCtx, util.magCtx, start,
                           util.resolutions, util.rng);

    // Extract measurements from result of propagation.
    KalmanFitterInputTrajectory traj;
    traj.startParameters = start;
    traj.sourceLinks = measurements.sourceLinks;

    trajectories.push_back(std::move(traj));
  }
  return trajectories;
}
}  // namespace ActsTests::AlignmentUtils

///
