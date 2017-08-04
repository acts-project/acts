// This file is part of the ACTS project.
//
// Copyright (C) 2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// STL include(s)
#include <memory>

// Boost include(s)
#define BOOST_TEST_MODULE Measurement Tests
#include <boost/test/included/unit_test.hpp>

// ATS include(s)
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/EventData/Measurement.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Fitter/KalmanUpdator.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Tools/LayerArrayCreator.hpp"
#include "ACTS/Tools/LayerCreator.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArrayXD.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

#include "KalmanFilterTestUtils.hpp"

/// use extrapolation to generate some measurements
std::vector<FitMeas_t>
generateDummyMeasurements(
    BoundParameters                               theTrackParameters,
    std::shared_ptr<IExtrapolationEngine>         theExtrapolationEngine,
    std::shared_ptr<const Acts::TrackingGeometry> geo)
{

  ExtrapolationCell<TrackParameters> exCell(theTrackParameters);
  exCell.addConfigurationMode(ExtrapolationMode::CollectSensitive);
  theExtrapolationEngine->extrapolate(exCell);

  std::vector<FitMeas_t> vMeasurements;
  vMeasurements.reserve(exCell.extrapolationSteps.size());

  // identifier
  long int id = 0;

  for (const auto& step : exCell.extrapolationSteps) {
    /// @TODO: accessing the local coordinates on the last extrapolation
    /// step seems to result in a segfault
    if (id >= exCell.extrapolationSteps.size() - 1) continue;
    const auto& tp = step.parameters;

    double            std1 = 0.01;
    double            std2 = 0.01;
    double            l1   = tp->get<eLOC_0>();
    double            l2   = tp->get<eLOC_1>();
    ActsSymMatrixD<2> cov;
    cov << std1 * std1, 0, 0, std2 * std2;

    vMeasurements.push_back(Meas_t<eLOC_0, eLOC_1>(
        tp->referenceSurface(), id, std::move(cov), l1, l2));
    ++id;
  }

  return vMeasurements;
}

namespace tt = boost::test_tools;
namespace Acts {
/// simple cylinder geometry for test purposes
std::shared_ptr<Acts::TrackingGeometry>
buildSimpleBarrel()
{
  // most placements are at the origin, without translations or rotations
  // create identity transformation to be used throughout
  std::shared_ptr<Acts::Transform3D> identityTransform
      = std::make_shared<Acts::Transform3D>();
  identityTransform->setIdentity();

  std::vector<Acts::LayerPtr> layerVector;
  // add barrel layers at hardcoded radii
  std::vector<double> layerRs{20, 40, 60, 80, 100};
  for (double r : layerRs) {
    Acts::Translation3D                translation{0., 0, 0};
    std::shared_ptr<Acts::Transform3D> transform
        = std::make_shared<Acts::Transform3D>(translation);
    std::shared_ptr<Acts::CylinderBounds> radialBounds
        = std::make_shared<Acts::CylinderBounds>(r, 1000.);
    const Acts::Surface* cylinderSurface
        = new Acts::CylinderSurface(identityTransform, r, 500);
    auto bp = std::make_unique<BinnedArrayXD<const Surface*>>(cylinderSurface);
    Acts::LayerPtr cylinderLayer = Acts::CylinderLayer::create(
        transform, radialBounds, std::move(bp), 0, nullptr, Acts::passive);
    layerVector.push_back(cylinderLayer);
  }
  // Module material - X0, L0, A, Z, Rho
  std::shared_ptr<Acts::Material> material
      = std::make_shared<Acts::Material>(95.7, 465.2, 28.03, 14., 2.32e-3);
  std::shared_ptr<Acts::CylinderVolumeBounds> cylinderVolumeBounds
      = std::make_shared<Acts::CylinderVolumeBounds>(120, 1000);
  Acts::LayerArrayCreator lc;
  auto la = lc.layerArray(layerVector, 0, 1020, Acts::arbitrary, Acts::binR);
  Acts::TrackingVolumePtr trackingVolume
      = Acts::TrackingVolume::create(identityTransform,
                                     cylinderVolumeBounds,
                                     nullptr,
                                     std::move(la),
                                     {},
                                     {},
                                     {},
                                     "MyVolume");
  std::const_pointer_cast<TrackingVolume>(trackingVolume)->sign(Acts::Global);
  std::shared_ptr<Acts::TrackingGeometry> trackingGeometry
      = std::make_shared<Acts::TrackingGeometry>(
          std::const_pointer_cast<TrackingVolume>(trackingVolume));
  return trackingGeometry;
}

namespace Test {

  /// test fit of a known track
  BOOST_AUTO_TEST_CASE(KalmanFilterFromPerfectInitial)
  {
    auto           geo   = buildSimpleBarrel();
    const Surface* pSurf = geo->getBeamline();
    double         x     = 0;
    double         y     = 0;
    double         z     = 0;
    double         px    = 100;
    double         py    = 0;
    double         pz    = 0;
    double         q     = 1;
    Vector3D       pos(x, y, z);
    Vector3D       mom(px, py, pz);

    // start covariance matrix
    auto startCov = std::make_unique<ActsSymMatrix<ParValue_t, NGlobalPars>>(
        ActsSymMatrix<ParValue_t, NGlobalPars>::Identity());
    (*startCov) = (*startCov) * 0.0001;

    auto startTP = std::make_unique<BoundParameters>(
        std::move(startCov), std::move(pos), std::move(mom), q, *pSurf);

    auto exEngine = initExtrapolator(geo);

    auto vMeasurements = generateDummyMeasurements(*startTP, exEngine, geo);

    KalmanFitter<MyExtrapolator,
                 CacheGenerator,
                 NoCalibration,
                 GainMatrixUpdator>
        KF;
    KF.m_oCacheGenerator = CacheGenerator();
    KF.m_oCalibrator     = NoCalibration();
    KF.m_oExtrapolator   = MyExtrapolator(exEngine);
    KF.m_oUpdator        = GainMatrixUpdator();
    auto track           = KF.fit(vMeasurements, std::move(startTP));

    int trackCounter = 0;
    for (const auto& p : track) {
      auto smoothedState = *p->getSmoothedState();
      auto filteredState = *p->getFilteredState();
      /// Test that position obtained by smoothed and filtered state are the
      /// same (they should be
      /// because the initial state describes the track perfectly)
      BOOST_TEST(smoothedState.position().norm()
                     == filteredState.position().norm(),
                 tt::tolerance(1e-7));
      ++trackCounter;
    }
  }
}  // end of namespace Test
}  // end of namespace Acts
