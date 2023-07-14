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
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#include <vector>

#include "FitterTestsCommon.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

//// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters(const ActsScalar x=0.0_m,
                                                const ActsScalar y=0.0_m,
                                                const ActsScalar z=0.0_m,
                                                const ActsScalar w=42_ns,
                                                const ActsScalar phi=0_degree,
                                                const ActsScalar theta=90_degree,
                                                const ActsScalar p=1_GeV,
                                                const ActsScalar q=1_e) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(x, y, z, w);
  return Acts::CurvilinearTrackParameters(mPos4, phi, theta, p, q, cov);
}

// Construct a straight-line propagator.
auto makeStraightPropagator(std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{std::move(geo)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::StraightLineStepper stepper;
  return Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>(
      stepper, std::move(navigator));
}

static std::vector<Acts::SourceLink> prepareSourceLinks(
    const std::vector<TestSourceLink>& sourceLinks) {
  std::vector<Acts::SourceLink> result;
  std::transform(sourceLinks.begin(), sourceLinks.end(),
                 std::back_inserter(result),
                 [](const auto& sl) { return Acts::SourceLink{sl}; });
  return result;
}

std::shared_ptr<const TrackingGeometry> makeToyDetector(
    const GeometryContext& tgContext, const size_t nSurfaces = 5) {
  if (nSurfaces < 1) {
    throw std::invalid_argument("At least 1 surfaces needs to be created.");
  }
  // Construct builder
  CuboidVolumeBuilder cvb;

  // Create configurations for surfaces
  std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
  for (unsigned int i = 1; i <= nSurfaces; i++) {
    // Position of the surfaces
    CuboidVolumeBuilder::SurfaceConfig cfg;
    cfg.position = {i * UnitConstants::m, 0, 0.};

    // Rotation of the surfaces
    double rotationAngle = M_PI * 0.5;
    Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3 yPos(0., 1., 0.);
    Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    cfg.rotation.col(0) = xPos;
    cfg.rotation.col(1) = yPos;
    cfg.rotation.col(2) = zPos;
    /// Shape of the surface
    // Boundaries of the surfaces
    cfg.rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));

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
    layerConfig.push_back(cfg);
  }

  for (auto& cfg : layerConfig) {
    cfg.surfaces = {};
  }

  // Inner Volume - Build volume configuration
  CuboidVolumeBuilder::VolumeConfig volumeConfig;
  volumeConfig.position = {2.5_m, 0., 0.};
  volumeConfig.length = {5_m, 1_m, 1_m};
  volumeConfig.layerCfg = layerConfig;
  volumeConfig.name = "Test volume";
  volumeConfig.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeBeryllium());

  volumeConfig.layers.clear();
  for (auto& lay : volumeConfig.layerCfg) {
    lay.active = true;
  }

  // Outer volume - Build TrackingGeometry configuration
  CuboidVolumeBuilder::Config config;
  config.position = {2.5_m, 0., 0.};
  config.length = {5_m, 1_m, 1_m};
  config.volumeCfg = {volumeConfig};

  cvb.setConfig(config);

  TrackingGeometryBuilder::Config tgbCfg;

  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });

  TrackingGeometryBuilder tgb(tgbCfg);

  std::unique_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);
  return detector;
}

struct Detector {
  // geometry
  std::shared_ptr<const TrackingGeometry> geometry;
};

BOOST_AUTO_TEST_SUITE(Gx2fTest)

// This test checks if the call to the fitter works and no errors occur in the
// framework, without fitting and updating any parameters
BOOST_AUTO_TEST_CASE(NoFit) {
  std::cout << "\n\n##### Start test case NoFit #####" << std::endl;

  // Context objects
  Acts::GeometryContext geoCtx;
  Acts::MagneticFieldContext magCtx;
  Acts::CalibrationContext calCtx;
  std::default_random_engine rng(42);

  Detector detector;
  const size_t nSurfaces = 5;
  detector.geometry = makeToyDetector(geoCtx, nSurfaces);

  auto parametersMeasurements = makeParameters();
  auto startParametersFit = makeParameters(0.1_m, 0.1_m, 0.1_m, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);

  MeasurementResolution resPixel = {MeasurementType::eLoc01, {25_um, 50_um}};
  MeasurementResolutionMap resolutions = {
      {Acts::GeometryIdentifier().setVolume(0), resPixel}};

  // propagator
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, parametersMeasurements,
                                         resolutions, rng);
  auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);

  using Gx2Fitter =
      Experimental::Gx2Fitter<SimPropagator, VectorMultiTrajectory>;
  Gx2Fitter Fitter(simPropagator);

  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&Test::testSourceLinkCalibrator<VectorMultiTrajectory>>();

  Experimental::Gx2FitterOptions gx2fOptions(
      geoCtx, magCtx, calCtx, extensions, PropagatorPlainOptions(), rSurface,
      false, false, FreeToBoundCorrection(false), 0);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};

  // Fit the track
  auto res = Fitter.fit(sourceLinks.begin(), sourceLinks.end(), startParametersFit,
                        gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), Acts::MultiTrajectoryTraits::kInvalid);
  BOOST_CHECK(!track.hasReferenceSurface());
  BOOST_CHECK_EQUAL(track.nMeasurements(), 0u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_EQUAL(track.parameters(), startParametersFit.parameters());
  BOOST_CHECK_EQUAL(track.covariance(), BoundMatrix::Identity());

  std::cout << "##### Finished test case NoFit #####" << std::endl;
}

BOOST_AUTO_TEST_CASE(Fit5Iterations) {
  std::cout << "\n\n##### Start test case Fit5Iterations #####" << std::endl;

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  Detector detector;
  const size_t nSurfaces = 5;
  detector.geometry = makeToyDetector(tgContext, nSurfaces);

  std::cout << "\n*** Go to propagator ***\n" << std::endl;

  auto parametersMeasurements = makeParameters();
  auto startParametersFit = makeParameters(10_mm, 10_mm, 10_mm, 42_ns, 10_degree, 80_degree, 1_GeV, 1_e);
  //  auto startParametersFit = parametersMeasurements;
  // Context objects
  Acts::GeometryContext geoCtx;
  Acts::MagneticFieldContext magCtx;
  // Acts::CalibrationContext calCtx;
  std::default_random_engine rng(42);

  MeasurementResolution resPixel = {MeasurementType::eLoc01, {25_um, 50_um}};
  // MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
  // MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
  MeasurementResolutionMap resolutions = {
      {Acts::GeometryIdentifier().setVolume(0), resPixel}};

  // simulation propagator
  using SimPropagator =
      Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
  SimPropagator simPropagator = makeStraightPropagator(detector.geometry);
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx,
                                         parametersMeasurements,
                                         resolutions, rng);
  auto sourceLinks = prepareSourceLinks(measurements.sourceLinks);
  std::cout << "sourceLinks.size() = " << sourceLinks.size() << std::endl;

  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nSurfaces);

  std::cout << "\n*** Start fitting ***\n" << std::endl;

  std::cout << "\n*** startParameter unsmeared: ***\n" << parametersMeasurements << std::endl;
  std::cout << "\n*** startParameter fit: ***\n" << startParametersFit << std::endl;
  const Surface* rSurface = &parametersMeasurements.referenceSurface();

  Navigator::Config cfg{detector.geometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Navigator rNavigator(cfg);
  // Configure propagation with deactivated B-field
  auto bField = std::make_shared<ConstantBField>(Vector3(0., 0., 0.));
  using RecoStepper = EigenStepper<>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  using Gx2Fitter =
      Experimental::Gx2Fitter<RecoPropagator, VectorMultiTrajectory>;

  Gx2Fitter Fitter(rPropagator);

  Experimental::Gx2FitterExtensions<VectorMultiTrajectory> extensions;
  extensions.calibrator
      .connect<&Test::testSourceLinkCalibrator<VectorMultiTrajectory>>();

  MagneticFieldContext mfContext = MagneticFieldContext();
  CalibrationContext calContext = CalibrationContext();

  Experimental::Gx2FitterOptions gx2fOptions(
      tgContext, mfContext, calContext, extensions, PropagatorPlainOptions(),
      rSurface, false, false, FreeToBoundCorrection(false), 5);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                  Acts::VectorMultiTrajectory{}};

  // Fit the track
  auto res = Fitter.fit(sourceLinks.begin(), sourceLinks.end(), startParametersFit,
                        gx2fOptions, tracks);

  BOOST_REQUIRE(res.ok());

  auto& track = *res;
  BOOST_CHECK_EQUAL(track.tipIndex(), Acts::MultiTrajectoryTraits::kInvalid);
  BOOST_CHECK(!track.hasReferenceSurface());
  BOOST_CHECK_EQUAL(track.nMeasurements(), 0u);
  BOOST_CHECK_EQUAL(track.nHoles(), 0u);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc0], -10.521387980023833, 1e-8);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundLoc1], -10.143918540958712, 1e-8);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundPhi], 1.4992306159443798e-05, 1e-8);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTheta], 1.5707952030421313, 1e-8);
  BOOST_CHECK_EQUAL(track.parameters()[eBoundQOverP], 1);
  BOOST_CHECK_CLOSE(track.parameters()[eBoundTime], 12591.283236000001, 1e-8);
  BOOST_CHECK_CLOSE(track.covariance().determinant(), 1.0382256071672664e-27, 1e-8);

  std::cout << "##### Finished test case Fit5Iterations #####" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
