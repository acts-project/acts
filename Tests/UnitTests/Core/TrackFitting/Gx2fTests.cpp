// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <vector>

#include "FitterTestsCommon.hpp"

using namespace Acts::UnitLiterals;
using namespace Acts::detail::Test;

Acts::Logging::Level logLevel = Acts::Logging::VERBOSE;
const auto gx2fLogger = Acts::getDefaultLogger("Gx2f", logLevel);

namespace Acts::Test {

//// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters(
    const ActsScalar x = 0.0_m, const ActsScalar y = 0.0_m,
    const ActsScalar z = 0.0_m, const ActsScalar w = 42_ns,
    const ActsScalar phi = 0_degree, const ActsScalar theta = 90_degree,
    const ActsScalar p = 2_GeV, const ActsScalar q = 1_e) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  const Acts::BoundSquareMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  const Acts::Vector4 mPos4(x, y, z, w);
  return Acts::CurvilinearTrackParameters(mPos4, phi, theta, q / p, cov,
                                          Acts::ParticleHypothesis::pion());
}

static std::vector<Acts::SourceLink> prepareSourceLinks(
    const std::vector<TestSourceLink>& sourceLinks) {
  std::vector<Acts::SourceLink> result;
  std::transform(sourceLinks.begin(), sourceLinks.end(),
                 std::back_inserter(result),
                 [](const auto& sl) { return Acts::SourceLink{sl}; });
  return result;
}

/// @brief Create a simple telescope detector
///
/// @param geoCtx
/// @param nSurfaces Number of surfaces
std::shared_ptr<const TrackingGeometry> makeToyDetector(
    const Acts::GeometryContext& geoCtx, const std::size_t nSurfaces = 5) {
  if (nSurfaces < 1) {
    throw std::invalid_argument("At least 1 surfaces needs to be created.");
  }

  // Define the dimensions of the square surfaces
  const double halfSizeSurface = 1_m;

  // Rotation of the surfaces
  const double rotationAngle = M_PI * 0.5;
  const Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
  const Vector3 yPos(0., 1., 0.);
  const Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));

  // Construct builder
  CuboidVolumeBuilder cvb;

  // Create configurations for surfaces
  std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
  for (unsigned int surfPos = 1; surfPos <= nSurfaces; surfPos++) {
    // Position of the surfaces
    CuboidVolumeBuilder::SurfaceConfig cfg;
    cfg.position = {surfPos * UnitConstants::m, 0., 0.};

    // Rotation of the surfaces
    cfg.rotation.col(0) = xPos;
    cfg.rotation.col(1) = yPos;
    cfg.rotation.col(2) = zPos;

    // Boundaries of the surfaces (shape)
    cfg.rBounds = std::make_shared<const RectangleBounds>(
        RectangleBounds(halfSizeSurface, halfSizeSurface));

    // Material of the surfaces
    MaterialSlab matProp(makeBeryllium(), 0.5_mm);
    cfg.surMat = std::make_shared<HomogeneousSurfaceMaterial>(matProp);

    // Thickness of the detector element
    cfg.thickness = 1_um;

    cfg.detElementConstructor =
        [](const Transform3& trans,
           const std::shared_ptr<const RectangleBounds>& bounds,
           double thickness) {
          return new DetectorElementStub(trans, bounds, thickness);
        };
    surfaceConfig.push_back(cfg);
  }

  // Build layer configurations
  std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig;
  for (auto& sCfg : surfaceConfig) {
    CuboidVolumeBuilder::LayerConfig cfg;
    cfg.surfaceCfg = {sCfg};
    cfg.active = true;
    cfg.envelopeX = {-0.1_mm, 0.1_mm};
    cfg.envelopeY = {-0.1_mm, 0.1_mm};
    cfg.envelopeZ = {-0.1_mm, 0.1_mm};
    layerConfig.push_back(cfg);
  }

  // Inner Volume - Build volume configuration
  CuboidVolumeBuilder::VolumeConfig volumeConfig;
  volumeConfig.length = {(nSurfaces + 1) * 1_m, 2 * halfSizeSurface,
                         2 * halfSizeSurface};
  volumeConfig.position = {volumeConfig.length.x() / 2, 0., 0.};
  volumeConfig.layerCfg = layerConfig;
  volumeConfig.name = "Test volume";

  // Outer volume - Build TrackingGeometry configuration
  CuboidVolumeBuilder::Config config;
  config.length = {(nSurfaces + 1) * 1_m, 2 * halfSizeSurface,
                   2 * halfSizeSurface};
  config.position = {volumeConfig.length.x() / 2, 0., 0.};
  config.volumeCfg = {volumeConfig};

  cvb.setConfig(config);

  TrackingGeometryBuilder::Config tgbCfg;

  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });

  TrackingGeometryBuilder tgb(tgbCfg);

  std::unique_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(geoCtx);
  return detector;
}

struct Detector {
  // geometry
  std::shared_ptr<const TrackingGeometry> geometry;
};

BOOST_AUTO_TEST_SUITE(Gx2fTest)
ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("Gx2fTests", logLevel))

// Context objects
const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;
const Acts::CalibrationContext calCtx;

// Measurement resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {25_um}};
const MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {50_um}};
const MeasurementResolutionMap resMapAllPixel = {
    {Acts::GeometryIdentifier().setVolume(0), resPixel}};

// This test checks if the call to the fitter works and no errors occur in the
// framework, without fitting and updating any parameters
BOOST_AUTO_TEST_CASE(NoFit) {
  ACTS_INFO("*** Test: NoFit -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 5;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  ACTS_DEBUG("Set up the fitter");
  // Reuse the SimPropagator, since we will not actually use it
  using Gx2Fitter =
      Experimental::Gx2Fitter<SimPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(simPropagator, gx2fLogger->clone());

  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 0, true, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), Acts::MultiTrajectoryTraits::kInvalid);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  BOOST_CHECK_EQUAL(track.chi2(), 0.);
  BOOST_CHECK_EQUAL(track.nDoF(), 0u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), 0u);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  BOOST_CHECK_EQUAL(track.parameters(), startParametersFit.parameters());
  BOOST_CHECK_EQUAL(track.covariance(), BoundMatrix::Identity());

  // Convergence
  BOOST_CHECK_EQUAL(
      (track.template component<
          std::size_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      0);

  ACTS_INFO("*** Test: NoFit -- Finish");
}

BOOST_AUTO_TEST_CASE(Fit5Iterations) {
  ACTS_INFO("*** Test: Fit5Iterations -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 5;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 5, true, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;

  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 8., 2.);
  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.2832360000, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1e-27, 4e0);

  // Convergence
  BOOST_CHECK_EQUAL(
      (track.template component<
          std::size_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      5);

  ACTS_INFO("*** Test: Fit5Iterations -- Finish");
}

BOOST_AUTO_TEST_CASE(MixedDetector) {
  ACTS_INFO("*** Test: MixedDetector -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 7;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  const MeasurementResolutionMap resMap = {
      {Acts::GeometryIdentifier().setVolume(2).setLayer(2), resPixel},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(4), resStrip0},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(6), resStrip1},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(8), resPixel},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(10), resStrip0},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(12), resStrip1},
      {Acts::GeometryIdentifier().setVolume(2).setLayer(14), resPixel},
  };

  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements = createMeasurements(
      simPropagator, geoCtx, magCtx, parametersMeasurements, resMap, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 5, true, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 8.5, 4.);
  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.2832360000, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 2e-28, 1e0);

  // Convergence
  BOOST_CHECK_EQUAL(
      (track.template component<
          std::size_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      5);

  ACTS_INFO("*** Test: MixedDetector -- Finish");
}

// This test checks if we can fit QOverP, when a magnetic field is introduced
BOOST_AUTO_TEST_CASE(FitWithBfield) {
  ACTS_INFO("*** Test: FitWithBfield -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 5;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  using SimStepper = EigenStepper<>;
  const auto simPropagator =
      makeConstantFieldPropagator<SimStepper>(detector.geometry, 0.3_T);

  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);

  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  // Reuse the SimPropagator, since it already uses the EigenStepper<>
  using SimPropagator = decltype(simPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<SimPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(simPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 5, false, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 7.5, 1.5);
  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  // TODO investigate further the reference values for eBoundPhi and det(cov)
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 8e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-4, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundQOverP], 0.5, 2e-1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.2832360000, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 8e-35, 4e0);

  // Convergence
  BOOST_CHECK_EQUAL(
      (track.template component<
          std::size_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      5);

  ACTS_INFO("*** Test: FitWithBfield -- Finish");
}

BOOST_AUTO_TEST_CASE(relChi2changeCutOff) {
  ACTS_INFO("*** Test: relChi2changeCutOff -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 5;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  // simulation propagator
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 500, true, 1e-5);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 8., 2.);
  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nSurfaces);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.2832360000, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1e-27, 4e0);

  // Convergence
  BOOST_CHECK_LT(
      (track.template component<
          std::size_t,
          hashString(Experimental::Gx2fConstants::gx2fnUpdateColumn)>()),
      10);

  ACTS_INFO("*** Test: relChi2changeCutOff -- Finish");
}

BOOST_AUTO_TEST_CASE(DidNotConverge) {
  ACTS_INFO("*** Test: DidNotConverge -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 5;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  // simulation propagator
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  // The relChi2changeCutOff = 0 prevents to stop the fitter after convergence,
  // therefore all updates will be done (even if the result does not change).
  // Since we didn't break due to convergence, we reach nUpdatesMax and
  // therefore fail the fit.
  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 6, true, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(!res.ok());
  BOOST_CHECK_EQUAL(
      res.error(),
      Acts::Experimental::GlobalChiSquareFitterError::DidNotConverge);

  ACTS_INFO("*** Test: DidNotConverge -- Finish");
}

BOOST_AUTO_TEST_CASE(NotEnoughMeasurements) {
  ACTS_INFO("*** Test: NotEnoughMeasurements -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  const std::size_t nSurfaces = 2;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  // simulation propagator
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  const auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() = " << sourceLinks.size());

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 6, true, 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(!res.ok());
  BOOST_CHECK_EQUAL(
      res.error(),
      Acts::Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements);

  ACTS_INFO("*** Test: NotEnoughMeasurements -- Finish");
}

BOOST_AUTO_TEST_CASE(FindHoles) {
  ACTS_INFO("*** Test: FindHoles -- Start");

  std::default_random_engine rng(42);

  ACTS_DEBUG("Create the detector");
  //  const std::size_t nSurfaces = 7;
  const std::size_t nSurfaces = 8;
  Detector detector;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  ACTS_DEBUG("Set the start parameters for measurement creation and fit");
  const auto parametersMeasurements = makeParameters();
  const auto startParametersFit = makeParameters(
      7_mm, 11_mm, 15_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  ACTS_DEBUG("Create the measurements");
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  const SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  const auto measurements =
      createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                         resMapAllPixel, rng);
  auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  ACTS_VERBOSE("sourceLinks.size() [before] = " << sourceLinks.size());

  // We remove the first measurement in the list. This does not create a hole.
  sourceLinks.erase(std::next(sourceLinks.begin(), 0));
  ACTS_VERBOSE(
      "sourceLinks.size() [after first erase] = " << sourceLinks.size());

  // We remove the last measurement in the list. This does not create a hole.
  sourceLinks.pop_back();
  ACTS_VERBOSE("sourceLinks.size() [after pop] = " << sourceLinks.size());

  // We remove the second to last measurement in the list. This effectivly
  // creates a hole on that surface.
  const std::size_t indexHole = sourceLinks.size() - 2;
  ACTS_VERBOSE("Remove measurement " << indexHole);
  sourceLinks.erase(std::next(sourceLinks.begin(), indexHole));
  ACTS_VERBOSE("sourceLinks.size() [after second-to-last erase]= "
               << sourceLinks.size());

  // We removed 3 measurements
  //  const std::size_t nMeasurements = nSurfaces - 2;
  const std::size_t nMeasurements = nSurfaces - 3;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  ACTS_DEBUG("Set up the fitter");
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  using RecoStepper = EigenStepper<>;
  const auto recoPropagator =
      makeConstantFieldPropagator<RecoStepper>(detector.geometry, 0_T);

  using RecoPropagator = decltype(recoPropagator);
  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;
  const Gx2Fitter fitter(recoPropagator, gx2fLogger->clone());

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&testSourceLinkCalibrator<VectorMultiTrajectory>>();
  TestSourceLink::SurfaceAccessor surfaceAccessor{*detector.geometry};
  extensions.surfaceAccessor
      .connect<&TestSourceLink::SurfaceAccessor::operator()>(&surfaceAccessor);

  const Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 20, true, 1e-5);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  ACTS_DEBUG("Fit the track");
  ACTS_VERBOSE("startParameter unsmeared:\n" << parametersMeasurements);
  ACTS_VERBOSE("startParameter fit:\n" << startParametersFit);
  const auto res = fitter.fit(sourceLinks.begin(), sourceLinks.end(),
                              startParametersFit, gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  const auto& track = *res;

  // -1, because the index starts at 0
  // -2, because the first and the last surface are not part of the track
  BOOST_CHECK_EQUAL(track.tipIndex(), nSurfaces - 1 - 2);
  BOOST_CHECK(track.hasReferenceSurface());

  // Track quantities
  CHECK_CLOSE_ABS(track.chi2(), 6.5, 2.);
  BOOST_CHECK_EQUAL(track.nDoF(), 10u);
  BOOST_CHECK_EQUAL(track.nHoles(), 1u);
  BOOST_CHECK_EQUAL(track.nMeasurements(), nMeasurements);
  BOOST_CHECK_EQUAL(track.nSharedHits(), 0u);
  BOOST_CHECK_EQUAL(track.nOutliers(), 0u);

  // Parameters
  // We need quite coarse checks here, since on different builds
  // the created measurements differ in the randomness
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -11., 7e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -15., 6e0);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1e-5, 1e3);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], M_PI / 2, 1e-3);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.2832360000, 1e-6);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 4.7e-28, 2e0);

  ACTS_INFO("*** Test: FindHoles -- Finish");
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace Acts::Test
