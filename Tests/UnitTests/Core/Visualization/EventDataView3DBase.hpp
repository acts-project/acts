// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Visualization/EventDataView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"

#include <cmath>
#include <fstream>
#include <optional>
#include <random>
#include <sstream>
#include <string>

using Acts::VectorHelpers::makeVector4;

namespace Acts {
namespace EventDataView3DTest {

using Covariance = BoundSymMatrix;
template <BoundIndices... params>
using MeasurementType =
    Measurement<Test::TestSourceLink, BoundIndices, params...>;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
///
/// @return an overall string including all written output
static inline std::string testBoundTrackParameters(IVisualization3D& helper) {
  std::stringstream ss;

  ViewConfig pcolor({20, 120, 20});
  ViewConfig scolor({235, 198, 52});

  auto gctx = GeometryContext();
  auto identity = Transform3::Identity();

  // rectangle and plane
  auto rectangle = std::make_shared<RectangleBounds>(15., 15.);
  auto plane = Surface::makeShared<PlaneSurface>(identity, rectangle);

  double momentumScale = 0.005;
  double localErrorScale = 10.;
  double directionErrorScale = 100.;

  // now create parameters on this surface
  // l_x, l_y, phi, theta, q/p (1/p), t
  std::array<double, 6> pars_array = {
      {-0.1234, 4.8765, 0.45, 0.128, 0.001, 21.}};

  BoundTrackParameters::ParametersVector pars =
      BoundTrackParameters::ParametersVector::Zero();
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  BoundSymMatrix cov = BoundSymMatrix::Zero();
  cov << 0.25, 0.0042, -0.00076, 6.156e-06, -2.11e-07, 0, 0.0042, 0.859,
      -0.000173, 0.000916, -4.017e-08, 0, -0.00076, -0.000173, 2.36e-04,
      -2.76e-07, 1.12e-08, 0, 6.15e-06, 0.000916, -2.76e-07, 8.84e-04,
      -2.85e-11, 0, -2.11 - 07, -4.017e-08, 1.123e-08, -2.85 - 11, 1.26e-10, 0,
      0, 0, 0, 0, 0, 1;

  EventDataView3D::drawBoundTrackParameters(
      helper, BoundTrackParameters(plane, pars, std::move(cov)), gctx,
      momentumScale, localErrorScale, directionErrorScale, pcolor, scolor);

  helper.write("EventData_BoundAtPlaneParameters");
  helper.write(ss);

  return ss.str();
}

static inline std::string testTrackState(IVisualization3D& helper) {
  std::stringstream ss;
  auto gctx = GeometryContext();

  double momentumScale = 0.1;
  double localErrorScale = 1.;
  double directionErrorScale = 1000;

  ViewConfig scolor({129, 173, 181});
  ViewConfig mcolor({255, 145, 48});
  mcolor.offset = -0.01;
  ViewConfig ppcolor({173, 133, 186});
  ppcolor.offset = -0.015;
  ViewConfig fpcolor({178, 200, 145});
  fpcolor.offset = -0.02;
  ViewConfig spcolor({20, 120, 20});
  spcolor.offset = -0.025;

  // Boundaries of the surfaces
  const auto rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(10, 10));
  // Material of the surfaces
  MaterialSlab matProp(Acts::Test::makeSilicon(), 0.5);
  const auto surfaceMaterial =
      std::make_shared<HomogeneousSurfaceMaterial>(matProp);
  // Construct the transform
  RotationMatrix3D rotation = RotationMatrix3D::Identity();
  double rotationAngle = M_PI / 2;
  Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3D yPos(0., 1., 0.);
  Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;
  Vector3D translation(0, 0, 0);
  // Vector3D translation(49, 0, 0);
  Transform3D trafo(Transform3D::Identity() * rotation);
  trafo.translation() = translation;
  // The surface
  auto surface = Surface::makeShared<PlaneSurface>(trafo, rBounds);
  // surface->assignSurfaceMaterial(surfaceMaterial);

  // The measurement
  const Vector2D loc(2.64678, 8.00062);
  ActsSymMatrixD<2> locCov;
  locCov << 0.0009, 0, 0, 0.0009;

  // The predicted parameters
  BoundVector predicted;
  BoundMatrix pCov;
  predicted << 2.58362, 7.88022, 0.010775, 1.5724, -0.25334, 140.272;
  pCov << 0.00280054, -1.16316e-07, 4.55413e-05, 2.0449e-09, -4.55897e-08,
      3.9782e-05, -1.16316e-07, 0.00280019, -2.34827e-09, -4.55322e-05,
      -1.02568e-12, 1.36071e-05, 4.55413e-05, -2.34827e-09, 8.38368e-07,
      3.73733e-11, -2.01704e-09, 6.4539e-07, 2.0449e-09, -4.55322e-05,
      3.73733e-11, 8.38124e-07, 4.1299e-14, -2.41901e-07, -4.55897e-08,
      -1.02568e-12, -2.01704e-09, 4.1299e-14, 6.4018e-08, -1.40676e-08,
      3.9782e-05, 1.36071e-05, 6.4539e-07, -2.41901e-07, -1.40676e-08,
      2.24689e+06;

  // The filtered parameters;
  BoundVector filtered;
  BoundMatrix fCov;
  filtered << 2.63142, 7.97133, 0.0115522, 1.57091, -0.253341, 140.273;
  fCov << 0.000681113, -6.88073e-09, 1.1076e-05, 1.49229e-10, -1.10878e-08,
      9.6754e-06, -6.88073e-09, 0.000681092, -2.22996e-10, -1.10748e-05,
      -5.98021e-13, 3.30997e-06, 1.1076e-05, -2.22996e-10, 2.77906e-07,
      9.258e-13, -1.45598e-09, 1.55809e-07, 1.49229e-10, -1.10748e-05,
      9.258e-13, 2.77834e-07, 3.62369e-14, -7.44673e-08, -1.10878e-08,
      -5.98021e-13, -1.45598e-09, 3.62369e-14, 6.40175e-08, -1.35774e-08,
      9.6754e-06, 3.30997e-06, 1.55809e-07, -7.44673e-08, -1.35774e-08,
      2.24689e+06;

  // The smootehd parameters;
  BoundVector smoothed;
  BoundMatrix sCov;
  smoothed << 2.63053, 7.98807, 0.0115394, 1.57052, -0.25334, 140.273;
  sCov << 0.000156495, 1.9583e-09, 3.59525e-07, 2.88542e-11, 4.57374e-08,
      2.25306e-06, 1.9583e-09, 0.000156456, -4.42585e-11, -3.58953e-07,
      -9.43293e-13, 2.98877e-07, 3.59525e-07, -4.42585e-11, 4.68873e-08,
      -6.69897e-13, 7.46465e-10, 4.1954e-09, 2.88542e-11, -3.58953e-07,
      -6.69897e-13, 4.68772e-08, -2.47294e-15, -1.04123e-08, 4.57374e-08,
      -9.43293e-13, 7.46465e-10, -2.47294e-15, 6.38984e-08, -1.27689e-08,
      2.25306e-06, 2.98877e-07, 4.1954e-09, -1.04123e-08, -1.27689e-08,
      2.24689e+06;

  EventDataView3D::drawTrackState(
      helper, gctx, *surface, loc, locCov,
      BoundTrackParameters(surface, predicted, std::move(pCov)),
      BoundTrackParameters(surface, filtered, std::move(fCov)),
      BoundTrackParameters(surface, smoothed, std::move(sCov)), momentumScale,
      localErrorScale, directionErrorScale, scolor, mcolor, ppcolor, fpcolor,
      spcolor);

  helper.write("EventData_TrackState");
  helper.write(ss);

  return ss.str();
}

static inline std::string testMultiTrajectory(IVisualization3D& helper) {
  std::stringstream ss;

  // Create a test context
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  CalibrationContext calContext = CalibrationContext();

  // Construct the rotation
  RotationMatrix3 rotation = RotationMatrix3::Identity();
  double rotationAngle = 90_degree;
  Vector3 xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3 yPos(0., 1., 0.);
  Vector3 zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Boundaries of the surfaces
  const auto rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(25_mm, 25_mm));

  // Material of the surfaces
  MaterialSlab matProp(Acts::Test::makeSilicon(), 0.5_mm);
  const auto surfaceMaterial =
      std::make_shared<HomogeneousSurfaceMaterial>(matProp);

  // Set translation vectors
  std::vector<Vector3> translations;
  translations.reserve(6);
  translations.push_back({-300_mm, 0., 0.});
  translations.push_back({-200_mm, 0., 0.});
  translations.push_back({-100_mm, 0., 0.});
  translations.push_back({100_mm, 0., 0.});
  translations.push_back({200_mm, 0., 0.});
  translations.push_back({300_mm, 0., 0.});

  // Construct layer configs
  std::vector<CuboidVolumeBuilder::LayerConfig> lConfs;
  lConfs.reserve(6);
  unsigned int i;
  for (i = 0; i < translations.size(); i++) {
    CuboidVolumeBuilder::SurfaceConfig sConf;
    sConf.position = translations[i];
    sConf.rotation = rotation;
    sConf.rBounds = rBounds;
    sConf.surMat = surfaceMaterial;
    // The thickness to construct the associated detector element
    sConf.thickness = 1_um;
    sConf.detElementConstructor =
        [](const Transform3& trans,
           std::shared_ptr<const RectangleBounds> bounds, double thickness) {
          return new Test::DetectorElementStub(trans, bounds, thickness);
        };
    CuboidVolumeBuilder::LayerConfig lConf;
    lConf.surfaceCfg = sConf;
    lConfs.push_back(lConf);
  }

  // Construct volume config
  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0., 0., 0.};
  vConf.length = {1.2_m, 1._m, 1._m};
  vConf.layerCfg = lConfs;
  vConf.name = "Tracker";

  // Construct volume builder config
  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {1.2_m, 1._m, 1._m};
  conf.volumeCfg = {vConf};  // one volume

  // Build detector
  std::cout << "Build the detector" << std::endl;
  CuboidVolumeBuilder cvb;
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  // Get the surfaces;
  std::vector<const Surface*> surfaces;
  surfaces.reserve(6);
  detector->visitSurfaces([&](const Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      std::cout << "surface " << surface->geometryId() << " placed at: ("
                << surface->center(tgContext).transpose() << " )" << std::endl;
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // Create measurements (assuming they are for a linear track parallel to
  // global x-axis)
  std::cout << "Creating measurements:" << std::endl;
  std::vector<Test::TestSourceLink> sourcelinks;
  sourcelinks.reserve(6);
  Vector2 lPosCenter{5_mm, 5_mm};
  Vector2 resolution{200_um, 150_um};
  SymMatrix2 cov2D = resolution.cwiseProduct(resolution).asDiagonal();
  for (const auto& surface : surfaces) {
    // 2D measurements
    Vector2 loc = lPosCenter;
    loc[0] += resolution[0] * gauss(generator);
    loc[1] += resolution[1] * gauss(generator);
    sourcelinks.emplace_back(eBoundLoc0, eBoundLoc1, loc, cov2D,
                             surface->geometryId());
  }

  // The KalmanFitter - we use the eigen stepper for covariance transport
  std::cout << "Construct KalmanFitter and perform fit" << std::endl;
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  // Set initial parameters for the particle track
  Covariance cov;
  cov << std::pow(100_um, 2), 0., 0., 0., 0., 0., 0., std::pow(100_um, 2), 0.,
      0., 0., 0., 0., 0., 0.0025, 0., 0., 0., 0., 0., 0., 0.0025, 0., 0., 0., 0.,
      0., 0., 0.01, 0., 0., 0., 0., 0., 0., 1.;
  Vector3D rPos(-350._mm, 100_um * gauss(generator), 100_um * gauss(generator));
  Vector3D rDir(1, 0.025 * gauss(generator), 0.025 * gauss(generator));
  CurvilinearTrackParameters rStart(makeVector4(rPos, 42_ns), rDir, 1_GeV, 1_e,
                                    cov);

  const Surface* rSurface = &rStart.referenceSurface();

  using Updater = GainMatrixUpdater;
  using Smoother = GainMatrixSmoother;
  using KalmanFitter = KalmanFitter<RecoPropagator, Updater, Smoother>;

  KalmanFitter kFitter(rPropagator);

  auto logger = getDefaultLogger("KalmanFilter", Logging::WARNING);
  KalmanFitterOptions<Test::TestSourceLinkCalibrator, VoidOutlierFinder>
      kfOptions(tgContext, mfContext, calContext,
                Test::TestSourceLinkCalibrator(), VoidOutlierFinder(),
                LoggerWrapper{*logger}, PropagatorPlainOptions(), rSurface);

  // Fit the track
  auto fitRes = kFitter.fit(sourcelinks, rStart, kfOptions);
  if (not fitRes.ok()) {
    std::cout << "Fit failed" << std::endl;
    return ss.str();
  }
  auto& fittedTrack = *fitRes;

  // Draw the track
  std::cout << "Draw the fitted track" << std::endl;
  double momentumScale = 10;
  double localErrorScale = 100.;
  double directionErrorScale = 100000;

  ViewConfig scolor({129, 173, 181});
  ViewConfig mcolor({255, 145, 48});
  mcolor.offset = -0.01;
  
  ViewConfig ppcolor({173, 133, 186});
  ppcolor.offset = -0.02;
  ViewConfig fpcolor({178, 200, 145});
  fpcolor.offset = -0.03;
  ViewConfig spcolor({20, 120, 20});
  spcolor.offset = -0.04;


  EventDataView3D::drawMultiTrajectory(
      helper, fittedTrack.fittedStates, fittedTrack.trackTip, tgContext,
      momentumScale, localErrorScale, directionErrorScale, scolor, mcolor,
      ppcolor, fpcolor, spcolor);

  helper.write("EventData_MultiTrajectory");
  helper.write(ss);

  return ss.str();
}

}  // namespace EventDataView3DTest
}  // namespace Acts
