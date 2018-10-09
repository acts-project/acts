// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// STL include(s)
#include <memory>

#define BOOST_TEST_MODULE KalmanFitter Tests
#include <boost/test/included/unit_test.hpp>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Fitter/KalmanActor.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Fitter/KalmanSequencer.hpp"
#include "Acts/Layers/PlaneLayer.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  // Propagator placeholder
  struct Propagator
  {
  };
  // Updator placeholder
  struct Updator
  {
  };
  // Calibrator placeholder
  struct Calibrator
  {
  };
  // Measurement placeholder
  struct MeasurementSet1
  {
    int ms = 1;
  };
  // Measurement placeholder
  struct MeasurementSet2
  {
    int ms = 2;
  };

  template <typename parameters_t>
  struct PropagatorState
  {
    // PropgatorState from parameters
    explicit PropagatorState(const parameters_t& pars) : stepping(pars) {}

    /// Fake a stepper state
    EigenStepper<ConstantBField>::State stepping;

    /// Fake a navigator state
    Navigator<KalmanSequencer>::State navigation;
  };

  // Shorthand
  using Jacobian   = ActsMatrixD<5, 5>;
  using Identifier = unsigned long int;
  template <ParID_t... params>
  using MeasurementType = Measurement<Identifier, params...>;
  template <ParID_t... params>
  using MeasuredState
      = MeasuredTrackState<Identifier, BoundParameters, Jacobian, params...>;
  using ParametricState
      = ParametricTrackState<Identifier, BoundParameters, Jacobian>;
  using VariantState = VariantTrackState<Identifier, BoundParameters, Jacobian>;
  using KalmanTrackStates = std::vector<VariantState>;

  // The plane surfaces
  PlaneSurface plane6(Vector3D(6., 0., 0.), Vector3D(1., 0., 0.));
  PlaneSurface plane7a(Vector3D(7., -0.1, 0.), Vector3D(1., 0., 0.));
  PlaneSurface plane7b(Vector3D(7., 0.1, 0.), Vector3D(1., 0., 0.));
  PlaneSurface plane8(Vector3D(8., 0., 0.), Vector3D(1., 0., 0.));
  PlaneSurface plane9(Vector3D(9., 0., 0.), Vector3D(1., 0., 0.));
  PlaneSurface plane10(Vector3D(10., 0., 0.), Vector3D(1., 0., 0.));

  auto planeLayer7 = PlaneLayer::create(nullptr, nullptr);
  auto planeLayer9 = PlaneLayer::create(nullptr, nullptr);

  /// Construct a KalmanActor
  BOOST_AUTO_TEST_CASE(KalmanActorConstruction)
  {

    plane7a.associateLayer(*planeLayer7);
    plane7b.associateLayer(*planeLayer7);
    plane9.associateLayer(*planeLayer9);

    // Construct the 1D measurement
    ActsSymMatrixD<1> cov1D7a;
    cov1D7a << 0.04;
    MeasurementType<ParDef::eLOC_0> m1D7a(plane7a, 0, std::move(cov1D7a), 0.02);

    ActsSymMatrixD<1> cov1D7b;
    cov1D7b << 0.04;
    MeasurementType<ParDef::eLOC_0> m1D7b(plane7b, 0, std::move(cov1D7b), 0.02);

    // Construct the 2D measurement
    ActsSymMatrixD<2> cov2D;
    cov2D << 0.04, 0., 0.09, 0.;
    MeasurementType<ParDef::eLOC_0, ParDef::eLOC_1> m2D(
        plane9, 0, std::move(cov2D), 0.02, 0.03);

    // The 1D track state from the measurement - on surface 7a
    VariantState mPlane7a = MeasuredState<ParDef::eLOC_0>(m1D7a);

    // The 1D track state from the measurement - on surface 7b
    VariantState mPlane7b = MeasuredState<ParDef::eLOC_0>(m1D7b);

    // The 2D track state from the measurement
    VariantState mPlane9 = MeasuredState<ParDef::eLOC_0, ParDef::eLOC_1>(m2D);

    using KalmanUpdator  = GainMatrixUpdator<BoundParameters, Jacobian>;
    using KalmanSmoother = GainMatrixSmoother<BoundParameters, Jacobian>;
    using KalmanActor
        = KalmanActor<KalmanTrackStates, KalmanUpdator, KalmanSmoother>;

    // Create the track states
    KalmanTrackStates tStates = {mPlane7a, mPlane7b, mPlane9};
    KalmanActor       kalmanActor;
    kalmanActor.trackStates = std::move(tStates);

    CurvilinearParameters cPars(
        nullptr, Vector3D(0., 0., 0.), Vector3D(0., 0., 0.), -1.);

    // The call objects
    PropagatorState<CurvilinearParameters> propState(cPars);
    KalmanActor::result_type               kalmanActorResult;
    kalmanActor(propState, kalmanActorResult);

    // Test that the fittedStates are now in the result
    BOOST_TEST(kalmanActorResult.fittedStates.size() == 3);

    // Test that the layer map of the PropagatorState is filled
    BOOST_TEST(propState.navigation.sequence.externalSurfaces.size() == 3);

    // Test that two measurements are on layer 7 & 9
    auto layer7ms = propState.navigation.sequence.externalSurfaces.count(
        planeLayer7.get());
    BOOST_TEST(layer7ms == 2);

    auto layer9ms = propState.navigation.sequence.externalSurfaces.count(
        planeLayer9.get());
    BOOST_TEST(layer9ms == 1);
  }

  /// Construct a KalmanFitter
  // BOOST_AUTO_TEST_CASE(KalmanFitterConstruction)
  //{
  //
  //
  //
  //  // Propgator, Updator and input measurements
  //  Propagator p;
  //  Updator u;
  //  Parameters ip;
  //  Surface s;
  //  MeasurementSet1 input1;
  //
  //  /// Testing the KF construction
  //  KalmanFitter<Propagator,Updator> kfitterM2M(p,u);
  //  auto output1 = kfitterM2M.fit(input1,ip,&s);
  //  BOOST_TEST(output1.ms == 1);
  //
  //  // Now create an input converter
  //
  //}

}  // namespace Test
}  // namespace Acts

/*#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Fitter/KalmanUpdator.hpp"
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArrayXD.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Volumes/CylinderVolumeBounds.hpp"

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
  size_t id = 0;
  for (const auto& step : exCell.extrapolationSteps) {
    /// @TODO: accessing the local coordinates on the last extrapolation
    /// step seems to result in a segfault
    if (id + 1 >= exCell.extrapolationSteps.size()) continue;
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
    auto           bp = std::make_unique<SurfaceArray>(cylinderSurface);
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
    double         x     = 0.;
    double         y     = 0.;
    double         z     = 0.;
    double         px    = 100.;
    double         py    = 0.;
    double         pz    = 0.;
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

      // Test that position obtained by smoothed and filtered state are
      // identical: they should be because the initial state describes
      // the track perfectly
      BOOST_TEST(smoothedState.position().norm()
                     == filteredState.position().norm(),
                 tt::tolerance(1e-7));

      ++trackCounter;
    }
  }
}  // namespace Test
}  // namespace Acts
*/