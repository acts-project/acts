// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
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
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/EventDataView.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"

#include <cmath>
#include <fstream>
#include <optional>
#include <random>
#include <sstream>
#include <string>

namespace Acts {
namespace EventDataViewTest {
using SourceLink = MinimalSourceLink;
using Covariance = BoundSymMatrix;

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
///
/// @return an overall string including all written output
static inline std::string testBoundParameters(IVisualization& helper) {
  std::stringstream ss;

  ViewConfig pcolor({20, 120, 20});
  ViewConfig scolor({235, 198, 52});

  auto gctx = GeometryContext();
  auto identity = std::make_shared<Transform3D>(Transform3D::Identity());

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

  BoundParameters::ParVector_t pars = BoundParameters::ParVector_t::Zero();
  pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
      pars_array[4], pars_array[5];

  BoundSymMatrix cov = BoundSymMatrix::Zero();
  cov << 0.25, 0.0042, -0.00076, 6.156e-06, -2.11e-07, 0, 0.0042, 0.859,
      -0.000173, 0.000916, -4.017e-08, 0, -0.00076, -0.000173, 2.36e-04,
      -2.76e-07, 1.12e-08, 0, 6.15e-06, 0.000916, -2.76e-07, 8.84e-04,
      -2.85e-11, 0, -2.11 - 07, -4.017e-08, 1.123e-08, -2.85 - 11, 1.26e-10, 0,
      0, 0, 0, 0, 0, 1;

  EventDataView::drawBoundParameters(
      helper, BoundParameters(gctx, std::move(cov), pars, plane), gctx,
      momentumScale, localErrorScale, directionErrorScale, pcolor, scolor);

  helper.write("EventData_BoundAtPlaneParameters");
  helper.write(ss);

  return ss.str();
}

static inline std::string testMultiTrajectory(IVisualization& helper) {
  std::stringstream ss;

  // Create a test context
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  CalibrationContext calContext = CalibrationContext();

  // Construct the rotation
  RotationMatrix3D rotation = RotationMatrix3D::Identity();
  double rotationAngle = 90_degree;
  Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3D yPos(0., 1., 0.);
  Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Boundaries of the surfaces
  const auto rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(0.1_m, 0.1_m));

  // Material of the surfaces
  MaterialProperties matProp(95.7, 465.2, 28.03, 14., 2.32e-3, 0.5_mm);
  const auto surfaceMaterial =
      std::make_shared<HomogeneousSurfaceMaterial>(matProp);

  // Set translation vectors
  std::vector<Vector3D> translations;
  translations.reserve(6);
  translations.push_back({-500_mm, 0., 0.});
  translations.push_back({-300_mm, 0., 0.});
  translations.push_back({-100_mm, 0., 0.});
  translations.push_back({100_mm, 0., 0.});
  translations.push_back({300_mm, 0., 0.});
  translations.push_back({500_mm, 0., 0.});

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
    sConf.thickness = 1._um;
    sConf.detElementConstructor =
        [](std::shared_ptr<const Transform3D> trans,
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
      std::cout << "surface " << surface->geoID() << " placed at: ("
                << surface->center(tgContext).transpose() << " )" << std::endl;
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // Create measurements (assuming they are for a linear track parallel to
  // global x-axis)
  std::cout << "Creating measurements:" << std::endl;
  std::vector<FittableMeasurement<SourceLink>> measurements;
  measurements.reserve(6);
  Vector2D lPosCenter{10_mm, 10_mm};
  std::array<double, 2> resolution = {30_um, 50_um};
  SymMatrix2D cov2D;
  cov2D << resolution[eLOC_0] * resolution[eLOC_0], 0., 0.,
      resolution[eLOC_1] * resolution[eLOC_1];
  for (const auto& surface : surfaces) {
    // 2D measurements
    double dx = resolution[eLOC_0] * gauss(generator);
    double dy = resolution[eLOC_1] * gauss(generator);
    MeasurementType<eLOC_0, eLOC_1> m01(surface->getSharedPtr(), {}, cov2D,
                                        lPosCenter[eLOC_0] + dx,
                                        lPosCenter[eLOC_1] + dy);
    measurements.push_back(std::move(m01));
  }

  // Make a vector of source links as input to the KF
  std::vector<SourceLink> sourcelinks;
  std::transform(measurements.begin(), measurements.end(),
                 std::back_inserter(sourcelinks),
                 [](const auto& m) { return SourceLink{&m}; });

  // The KalmanFitter - we use the eigen stepper for covariance transport
  std::cout << "Construct KalmanFitter and perform fit" << std::endl;
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  // Set initial parameters for the particle track
  Covariance cov;
  cov << std::pow(100_um, 2), 0., 0., 0., 0., 0., 0., std::pow(100_um, 2), 0.,
      0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0.,
      0., 0., 0.01, 0., 0., 0., 0., 0., 0., 1.;

  Vector3D rPos(-1_m, 100_um * gauss(generator), 100_um * gauss(generator));
  Vector3D rMom(1_GeV, 0.025_GeV * gauss(generator),
                0.025_GeV * gauss(generator));

  SingleCurvilinearTrackParameters<ChargedPolicy> rStart(cov, rPos, rMom, 1.,
                                                         42.);

  const Surface* rSurface = &rStart.referenceSurface();

  using Updater = GainMatrixUpdater;
  using Smoother = GainMatrixSmoother;
  using KalmanFitter = KalmanFitter<RecoPropagator, Updater, Smoother>;

  KalmanFitter kFitter(rPropagator,
                       getDefaultLogger("KalmanFilter", Logging::WARNING));

  KalmanFitterOptions<VoidOutlierFinder> kfOptions(
      tgContext, mfContext, calContext, VoidOutlierFinder(), rSurface);

  // Fit the track
  auto fitRes = kFitter.fit(sourcelinks, rStart, kfOptions);
  if (not fitRes.ok()) {
    std::cout << "Fit failed" << std::endl;
    return ss.str();
  }
  auto& fittedTrack = *fitRes;

  // Draw the track
  std::cout << "Draw the fitted track" << std::endl;
  double momentumScale = 15;
  double localErrorScale = 100.;
  double directionErrorScale = 500000;

  ViewConfig scolor({235, 198, 52});
  ViewConfig mcolor({255, 145, 48});
  mcolor.offset = -0.01;
  ViewConfig ppcolor({138, 214, 255});
  ppcolor.offset = -0.02;
  ViewConfig fpcolor({92, 149, 255});
  fpcolor.offset = -0.03;
  ViewConfig spcolor({20, 120, 20});
  spcolor.offset = -0.04;

  EventDataView::drawMultiTrajectory(
      helper, fittedTrack.fittedStates, fittedTrack.trackTip, tgContext,
      momentumScale, localErrorScale, directionErrorScale, scolor, mcolor,
      ppcolor, fpcolor, spcolor);

  helper.write("EventData_MultiTrajectory");
  helper.write(ss);

  return ss.str();
}

}  // namespace EventDataViewTest
}  // namespace Acts
