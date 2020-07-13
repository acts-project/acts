// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <fstream>

#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Propagator/DebugOutputActor.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

/// @brief Simplified propagator state
template <typename stepper_state_t>
struct PropState {
  /// @brief Constructor
  PropState(stepper_state_t sState) : stepping(sState) {}
  /// State of the eigen stepper
  stepper_state_t stepping;
  /// Propagator options which only carry the relevant components
  struct {
    double mass = 42.;
    double tolerance = 1e-4;
    double stepSizeCutOff = 0.;
    unsigned int maxRungeKuttaStepTrials = 10000;
  } options;
};

/// @brief Aborter for the case that a particle leaves the detector or reaches
/// a custom made threshold.
///
struct EndOfWorld {
  /// Maximum value in x-direction of the detector
  double maxX = 1_m;

  /// @brief Constructor
  EndOfWorld() = default;

  /// @brief Main call operator for the abort operation
  ///
  /// @tparam propagator_state_t State of the propagator
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  /// @return Boolean statement if the particle is still in the detector
  template <typename propagator_state_t, typename stepper_t>
  bool operator()(propagator_state_t& state, const stepper_t& stepper) const {
    const double tolerance = state.options.targetTolerance;
    if (maxX - std::abs(stepper.position(state.stepping).x()) <= tolerance ||
        std::abs(stepper.position(state.stepping).y()) >= 0.5_m ||
        std::abs(stepper.position(state.stepping).z()) >= 0.5_m)
      return true;
    return false;
  }
};

///
/// @brief Data collector while propagation
///
struct StepCollector {
  ///
  /// @brief Data container for result analysis
  ///
  struct this_result {
    // Position of the propagator after each step
    std::vector<Vector3D> position;
    // Momentum of the propagator after each step
    std::vector<Vector3D> momentum;
  };

  using result_type = this_result;

  /// @brief Main call operator for the action list. It stores the data for
  /// analysis afterwards
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [in] stepper Stepper of the propagation
  /// @param [out] result Struct which is filled with the data
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    result.position.push_back(stepper.position(state.stepping));
    result.momentum.push_back(stepper.momentum(state.stepping) *
                              stepper.direction(state.stepping));
  }
};

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(eigen_stepper_state_test) {
  // Set up some variables
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  ConstantBField bField(Vector3D(1., 2.5, 33.33));

  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  double time = 7.;
  double charge = -1.;

  // Test charged parameters without covariance matrix
  CurvilinearParameters cp(std::nullopt, pos, mom, charge, time);
  EigenStepper<ConstantBField>::State esState(tgContext, mfContext, cp, ndir,
                                              stepSize, tolerance);

  // Test the result & compare with the input/test for reasonable members
  BOOST_TEST(esState.jacToGlobal == BoundToFreeMatrix::Zero());
  BOOST_TEST(esState.jacTransport == FreeMatrix::Identity());
  BOOST_TEST(esState.derivative == FreeVector::Zero());
  BOOST_TEST(!esState.covTransport);
  BOOST_TEST(esState.cov == Covariance::Zero());
  BOOST_TEST(esState.pos == pos);
  BOOST_TEST(esState.dir == mom.normalized());
  BOOST_TEST(esState.p == mom.norm());
  BOOST_TEST(esState.q == charge);
  BOOST_TEST(esState.t == time);
  BOOST_TEST(esState.navDir == ndir);
  BOOST_TEST(esState.pathAccumulated == 0.);
  BOOST_TEST(esState.stepSize == ndir * stepSize);
  BOOST_TEST(esState.previousStepSize == 0.);
  BOOST_TEST(esState.tolerance == tolerance);

  // Test without charge and covariance matrix
  NeutralCurvilinearTrackParameters ncp(std::nullopt, pos, mom, time);
  esState = EigenStepper<ConstantBField>::State(tgContext, mfContext, ncp, ndir,
                                                stepSize, tolerance);
  BOOST_TEST(esState.q == 0.);

  // Test with covariance matrix
  Covariance cov = 8. * Covariance::Identity();
  ncp = NeutralCurvilinearTrackParameters(cov, pos, mom, time);
  esState = EigenStepper<ConstantBField>::State(tgContext, mfContext, ncp, ndir,
                                                stepSize, tolerance);
  BOOST_TEST(esState.jacToGlobal != BoundToFreeMatrix::Zero());
  BOOST_TEST(esState.covTransport);
  BOOST_TEST(esState.cov == cov);
}

/// These tests are aiming to test the functions of the EigenStepper
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(eigen_stepper_test) {
  // Set up some variables for the state
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  ConstantBField bField(Vector3D(1., 2.5, 33.33));

  // Construct the parameters
  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  double time = 7.;
  double charge = -1.;
  Covariance cov = 8. * Covariance::Identity();
  CurvilinearParameters cp(cov, pos, mom, charge, time);

  // Build the state and the stepper
  EigenStepper<ConstantBField>::State esState(tgContext, mfContext, cp, ndir,
                                              stepSize, tolerance);
  EigenStepper<ConstantBField> es(bField);

  // Test the getters
  BOOST_TEST(es.position(esState) == esState.pos);
  BOOST_TEST(es.direction(esState) == esState.dir);
  BOOST_TEST(es.momentum(esState) == esState.p);
  BOOST_TEST(es.charge(esState) == esState.q);
  BOOST_TEST(es.time(esState) == esState.t);
  //~ BOOST_TEST(es.overstepLimit(esState) == tolerance);
  BOOST_TEST(es.getField(esState, pos) == bField.getField(pos));

  // Step size modifies
  const std::string originalStepSize = esState.stepSize.toString();

  es.setStepSize(esState, 1337.);
  BOOST_TEST(esState.previousStepSize == ndir * stepSize);
  BOOST_TEST(esState.stepSize == 1337.);

  es.releaseStepSize(esState);
  BOOST_TEST(esState.stepSize == -123.);
  BOOST_TEST(es.outputStepSize(esState) == originalStepSize);

  // Test the curvilinear state construction
  auto curvState = es.curvilinearState(esState);
  auto curvPars = std::get<0>(curvState);
  CHECK_CLOSE_ABS(curvPars.position(), cp.position(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.momentum(), cp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), 1e-6);
  CHECK_CLOSE_ABS(curvPars.time(), cp.time(), 1e-6);
  BOOST_TEST(curvPars.covariance().has_value());
  BOOST_TEST(*curvPars.covariance() != cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(curvState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(curvState), 0., 1e-6);

  // Test the update method
  Vector3D newPos(2., 4., 8.);
  Vector3D newMom(3., 9., 27.);
  double newTime(321.);
  es.update(esState, newPos, newMom.normalized(), newMom.norm(), newTime);
  BOOST_TEST(esState.pos == newPos);
  BOOST_TEST(esState.dir == newMom.normalized());
  BOOST_TEST(esState.p == newMom.norm());
  BOOST_TEST(esState.q == charge);
  BOOST_TEST(esState.t == newTime);

  // The covariance transport
  esState.cov = cov;
  es.covarianceTransport(esState);
  BOOST_TEST(esState.cov != cov);
  BOOST_TEST(esState.jacToGlobal != BoundToFreeMatrix::Zero());
  BOOST_TEST(esState.jacTransport == FreeMatrix::Identity());
  BOOST_TEST(esState.derivative == FreeVector::Zero());

  // Perform a step without and with covariance transport
  esState.cov = cov;
  PropState ps(esState);

  ps.stepping.covTransport = false;
  double h = es.step(ps).value();
  BOOST_TEST(ps.stepping.stepSize == h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_TEST(ps.stepping.pos.norm() != newPos.norm());
  BOOST_TEST(ps.stepping.dir != newMom.normalized());
  BOOST_TEST(ps.stepping.q == charge);
  BOOST_TEST(ps.stepping.t < newTime);
  BOOST_TEST(ps.stepping.derivative == FreeVector::Zero());
  BOOST_TEST(ps.stepping.jacTransport == FreeMatrix::Identity());

  ps.stepping.covTransport = true;
  double h2 = es.step(ps).value();
  BOOST_TEST(h2 == h);
  CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
  BOOST_TEST(ps.stepping.pos.norm() != newPos.norm());
  BOOST_TEST(ps.stepping.dir != newMom.normalized());
  BOOST_TEST(ps.stepping.q == charge);
  BOOST_TEST(ps.stepping.t < newTime);
  BOOST_TEST(ps.stepping.derivative != FreeVector::Zero());
  BOOST_TEST(ps.stepping.jacTransport != FreeMatrix::Identity());

  /// Repeat with surface related methods
  auto plane = Surface::makeShared<PlaneSurface>(pos, mom.normalized());
  BoundParameters bp(tgContext, cov, pos, mom, charge, time, plane);
  esState = EigenStepper<ConstantBField>::State(tgContext, mfContext, cp, ndir,
                                                stepSize, tolerance);

  // Test the intersection in the context of a surface
  auto targetSurface = Surface::makeShared<PlaneSurface>(
      pos + ndir * 2. * mom.normalized(), mom.normalized());
  es.updateSurfaceStatus(esState, *targetSurface, BoundaryCheck(false));
  BOOST_TEST(esState.stepSize.value(ConstrainedStep::actor), ndir * 2.);

  // Test the step size modification in the context of a surface
  es.updateStepSize(
      esState,
      targetSurface->intersect(esState.geoContext, esState.pos,
                               esState.navDir * esState.dir, false),
      false);
  BOOST_TEST(esState.stepSize == 2.);
  esState.stepSize = ndir * stepSize;
  es.updateStepSize(
      esState,
      targetSurface->intersect(esState.geoContext, esState.pos,
                               esState.navDir * esState.dir, false),
      true);
  BOOST_TEST(esState.stepSize == 2.);

  // Test the bound state construction
  auto boundState = es.boundState(esState, *plane);
  auto boundPars = std::get<0>(boundState);
  CHECK_CLOSE_ABS(boundPars.position(), bp.position(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.momentum(), bp.momentum(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.charge(), bp.charge(), 1e-6);
  CHECK_CLOSE_ABS(boundPars.time(), bp.time(), 1e-6);
  BOOST_TEST(boundPars.covariance().has_value());
  BOOST_TEST(*boundPars.covariance() != cov);
  CHECK_CLOSE_COVARIANCE(std::get<1>(boundState),
                         BoundMatrix(BoundMatrix::Identity()), 1e-6);
  CHECK_CLOSE_ABS(std::get<2>(boundState), 0., 1e-6);

  // Update in context of a surface
  BoundParameters bpTarget(tgContext, 2. * cov, 2. * pos, 2. * mom,
                           -1. * charge, 2. * time, targetSurface);
  Vector3D dir = bpTarget.momentum().normalized();
  FreeVector freeParams;
  freeParams[eFreePos0] = bpTarget.position()[eX];
  freeParams[eFreePos1] = bpTarget.position()[eY];
  freeParams[eFreePos2] = bpTarget.position()[eZ];
  freeParams[eFreeTime] = bpTarget.time();
  freeParams[eFreeDir0] = dir[eMom0];
  freeParams[eFreeDir1] = dir[eMom1];
  freeParams[eFreeDir2] = dir[eMom2];
  freeParams[eFreeQOverP] = bpTarget.charge() / bpTarget.momentum().norm();

  es.update(esState, freeParams, *bpTarget.covariance());
  BOOST_TEST(esState.pos == 2. * pos);
  CHECK_CLOSE_ABS(esState.dir, mom.normalized(), 1e-6);
  BOOST_TEST(esState.p == 2. * mom.norm());
  BOOST_TEST(esState.q == 1. * charge);
  BOOST_TEST(esState.t == 2. * time);
  CHECK_CLOSE_COVARIANCE(esState.cov, Covariance(2. * cov), 1e-6);

  // Transport the covariance in the context of a surface
  es.covarianceTransport(esState, *plane);
  BOOST_TEST(esState.cov != cov);
  BOOST_TEST(esState.jacToGlobal != BoundToFreeMatrix::Zero());
  BOOST_TEST(esState.jacTransport == FreeMatrix::Identity());
  BOOST_TEST(esState.derivative == FreeVector::Zero());

  // Test a case where no step size adjustment is required
  ps.options.tolerance = 2. * 4.4258e+09;
  double h0 = esState.stepSize;
  es.step(ps);
  CHECK_CLOSE_ABS(h0, esState.stepSize, 1e-6);

  // Produce some errors
  NullBField nBfield;
  EigenStepper<NullBField> nes(nBfield);
  EigenStepper<NullBField>::State nesState(tgContext, mfContext, cp, ndir,
                                           stepSize, tolerance);
  PropState nps(nesState);
  // Test that we can reach the minimum step size
  nps.options.tolerance = 1e-21;
  nps.options.stepSizeCutOff = 1e20;
  auto res = nes.step(nps);
  BOOST_TEST(!res.ok());
  BOOST_TEST(res.error() == EigenStepperError::StepSizeStalled);

  // Test that the number of trials exceeds
  nps.options.stepSizeCutOff = 0.;
  nps.options.maxRungeKuttaStepTrials = 0.;
  res = nes.step(nps);
  BOOST_TEST(!res.ok());
  BOOST_TEST(res.error() == EigenStepperError::StepSizeAdjustmentFailed);
}

/// @brief This function tests the EigenStepper with the DefaultExtension and
/// the DenseEnvironmentExtension. The focus of this tests lies in the
/// choosing of the right extension for the individual use case. This is
/// performed with three different detectors:
/// a) Pure vaccuum -> DefaultExtension needs to act
/// b) Pure Be -> DenseEnvironmentExtension needs to act
/// c) Vacuum - Be - Vacuum -> Both should act and switch during the
/// propagation

// Test case a). The DenseEnvironmentExtension should state that it is not
// valid in this case.
BOOST_AUTO_TEST_CASE(step_extension_vacuum_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0.5_m, 0., 0.};
  vConf.length = {1_m, 1_m, 1_m};
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg.push_back(vConf);
  conf.position = {0.5_m, 0., 0.};
  conf.length = {1_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> vacuum =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviVac(vacuum);
  naviVac.resolvePassive = true;
  naviVac.resolveMaterial = true;
  naviVac.resolveSensitive = true;

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  Vector3D startParams(0., 0., 0.), startMom(1_GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(cov, startParams,
                                                       startMom, 1., 0.);

  // Create action list for surface collection
  ActionList<StepCollector> aList;
  AbortList<EndOfWorld> abortList;

  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                AbortList<EndOfWorld>>
      propOpts(tgContext, mfContext);
  propOpts.actionList = aList;
  propOpts.abortList = abortList;
  propOpts.maxSteps = 100;
  propOpts.maxStepSize = 1.5_m;

  // Build stepper and propagator
  ConstantBField bField(Vector3D(0., 0., 0.));
  EigenStepper<
      ConstantBField,
      StepperExtensionList<DefaultExtension, DenseEnvironmentExtension>,
      detail::HighestValidAuctioneer>
      es(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DefaultExtension,
                                               DenseEnvironmentExtension>,
                          detail::HighestValidAuctioneer>,
             Navigator>
      prop(es, naviVac);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Check that the propagation happend without interactions
  for (const auto& pos : stepResult.position) {
    CHECK_SMALL(pos.y(), 1_um);
    CHECK_SMALL(pos.z(), 1_um);
    if (pos == stepResult.position.back())
      CHECK_CLOSE_ABS(pos.x(), 1_m, 1_um);
  }
  for (const auto& mom : stepResult.momentum) {
    CHECK_CLOSE_ABS(mom, startMom, 1_keV);
  }

  // Rebuild and check the choice of extension
  ActionList<StepCollector> aListDef;

  // Set options for propagator
  PropagatorOptions<ActionList<StepCollector>, AbortList<EndOfWorld>>
      propOptsDef(tgContext, mfContext);
  propOptsDef.actionList = aListDef;
  propOptsDef.abortList = abortList;
  propOptsDef.maxSteps = 100;
  propOptsDef.maxStepSize = 1.5_m;

  EigenStepper<ConstantBField, StepperExtensionList<DefaultExtension>> esDef(
      bField);
  Propagator<
      EigenStepper<ConstantBField, StepperExtensionList<DefaultExtension>>,
      Navigator>
      propDef(esDef, naviVac);

  // Launch and collect results
  const auto& resultDef = propDef.propagate(sbtp, propOptsDef).value();
  const StepCollector::this_result& stepResultDef =
      resultDef.get<typename StepCollector::result_type>();

  // Check that the right extension was chosen
  // If chosen correctly, the number of elements should be identical
  BOOST_TEST(stepResult.position.size() == stepResultDef.position.size());
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.position[i], stepResultDef.position[i], 1_um);
  }
  BOOST_TEST(stepResult.momentum.size() == stepResultDef.momentum.size());
  for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.momentum[i], stepResultDef.momentum[i], 1_keV);
  }
}
// Test case b). The DefaultExtension should state that it is invalid here.
BOOST_AUTO_TEST_CASE(step_extension_material_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0.5_m, 0., 0.};
  vConf.length = {1_m, 1_m, 1_m};
  vConf.volumeMaterial = std::make_shared<const HomogeneousVolumeMaterial>(
      Material(352.8, 394.133, 9.012, 4., 1.848e-3));
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg.push_back(vConf);
  conf.position = {0.5_m, 0., 0.};
  conf.length = {1_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> material =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviMat(material);
  naviMat.resolvePassive = true;
  naviMat.resolveMaterial = true;
  naviMat.resolveSensitive = true;

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  Vector3D startParams(0., 0., 0.), startMom(5_GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(cov, startParams,
                                                       startMom, 1., 0.);

  // Create action list for surface collection
  ActionList<StepCollector> aList;
  AbortList<EndOfWorld> abortList;

  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                AbortList<EndOfWorld>>
      propOpts(tgContext, mfContext);
  propOpts.actionList = aList;
  propOpts.abortList = abortList;
  propOpts.maxSteps = 100;
  propOpts.maxStepSize = 1.5_m;
  propOpts.debug = true;

  // Build stepper and propagator
  ConstantBField bField(Vector3D(0., 0., 0.));
  EigenStepper<
      ConstantBField,
      StepperExtensionList<DefaultExtension, DenseEnvironmentExtension>,
      detail::HighestValidAuctioneer>
      es(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DefaultExtension,
                                               DenseEnvironmentExtension>,
                          detail::HighestValidAuctioneer>,
             Navigator>
      prop(es, naviMat);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Check that there occured interaction
  for (const auto& pos : stepResult.position) {
    CHECK_SMALL(pos.y(), 1_um);
    CHECK_SMALL(pos.z(), 1_um);
    if (pos == stepResult.position.front()) {
      CHECK_SMALL(pos.x(), 1_um);
    } else {
      BOOST_CHECK_GT(std::abs(pos.x()), 1_um);
    }
  }
  for (const auto& mom : stepResult.momentum) {
    CHECK_SMALL(mom.y(), 1_keV);
    CHECK_SMALL(mom.z(), 1_keV);
    if (mom == stepResult.momentum.front()) {
      CHECK_CLOSE_ABS(mom.x(), 5_GeV, 1_keV);
    } else {
      BOOST_CHECK_LT(mom.x(), 5_GeV);
    }
  }

  // Rebuild and check the choice of extension
  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                AbortList<EndOfWorld>>
      propOptsDense(tgContext, mfContext);
  propOptsDense.actionList = aList;
  propOptsDense.abortList = abortList;
  propOptsDense.maxSteps = 100;
  propOptsDense.maxStepSize = 1.5_m;
  propOptsDense.debug = true;

  // Build stepper and propagator
  EigenStepper<ConstantBField, StepperExtensionList<DenseEnvironmentExtension>>
      esDense(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DenseEnvironmentExtension>>,
             Navigator>
      propDense(esDense, naviMat);

  // Launch and collect results
  const auto& resultDense = propDense.propagate(sbtp, propOptsDense).value();
  const StepCollector::this_result& stepResultDense =
      resultDense.get<typename StepCollector::result_type>();

  // Check that the right extension was chosen
  // If chosen correctly, the number of elements should be identical
  BOOST_TEST(stepResult.position.size() == stepResultDense.position.size());
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.position[i], stepResultDense.position[i], 1_um);
  }
  BOOST_TEST(stepResult.momentum.size() == stepResultDense.momentum.size());
  for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
    CHECK_CLOSE_ABS(stepResult.momentum[i], stepResultDense.momentum[i], 1_keV);
  }

  ////////////////////////////////////////////////////////////////////

  // Re-launch the configuration with magnetic field
  bField.setField(0., 1_T, 0.);
  EigenStepper<
      ConstantBField,
      StepperExtensionList<DefaultExtension, DenseEnvironmentExtension>,
      detail::HighestValidAuctioneer>
      esB(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DefaultExtension,
                                               DenseEnvironmentExtension>,
                          detail::HighestValidAuctioneer>,
             Navigator>
      propB(esB, naviMat);

  const auto& resultB = propB.propagate(sbtp, propOptsDense).value();
  const StepCollector::this_result& stepResultB =
      resultB.get<typename StepCollector::result_type>();

  // Check that there occured interaction
  for (const auto& pos : stepResultB.position) {
    if (pos == stepResultB.position.front()) {
      CHECK_SMALL(pos, 1_um);
    } else {
      BOOST_CHECK_GT(std::abs(pos.x()), 1_um);
      CHECK_SMALL(pos.y(), 1_um);
      BOOST_CHECK_GT(std::abs(pos.z()), 1_um);
    }
  }
  for (const auto& mom : stepResultB.momentum) {
    if (mom == stepResultB.momentum.front()) {
      CHECK_CLOSE_ABS(mom, startMom, 1_keV);
    } else {
      BOOST_CHECK_NE(mom.x(), 5_GeV);
      CHECK_SMALL(mom.y(), 1_keV);
      BOOST_CHECK_NE(mom.z(), 0.);
    }
  }
}
// Test case c). Both should be involved in their part of the detector
BOOST_AUTO_TEST_CASE(step_extension_vacmatvac_test) {
  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::VolumeConfig vConfVac1;
  vConfVac1.position = {0.5_m, 0., 0.};
  vConfVac1.length = {1_m, 1_m, 1_m};
  vConfVac1.name = "First vacuum volume";
  CuboidVolumeBuilder::VolumeConfig vConfMat;
  vConfMat.position = {1.5_m, 0., 0.};
  vConfMat.length = {1_m, 1_m, 1_m};
  vConfMat.volumeMaterial = std::make_shared<const HomogeneousVolumeMaterial>(
      Material(352.8, 394.133, 9.012, 4., 1.848e-3));
  vConfMat.name = "Material volume";
  CuboidVolumeBuilder::VolumeConfig vConfVac2;
  vConfVac2.position = {2.5_m, 0., 0.};
  vConfVac2.length = {1_m, 1_m, 1_m};
  vConfVac2.name = "Second vacuum volume";
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg = {vConfVac1, vConfMat, vConfVac2};
  conf.position = {1.5_m, 0., 0.};
  conf.length = {3_m, 1_m, 1_m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> det = tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviDet(det);
  naviDet.resolvePassive = true;
  naviDet.resolveMaterial = true;
  naviDet.resolveSensitive = true;

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  Vector3D startParams(0., 0., 0.), startMom(5_GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(cov, startParams,
                                                       startMom, 1., 0.);

  // Create action list for surface collection
  AbortList<EndOfWorld> abortList;
  abortList.get<EndOfWorld>().maxX = 3_m;

  using DebugOutput = Acts::DebugOutputActor;

  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector, DebugOutput>,
                                AbortList<EndOfWorld>>
      propOpts(tgContext, mfContext);
  propOpts.abortList = abortList;
  propOpts.maxSteps = 100;
  propOpts.maxStepSize = 1.5_m;

  // Build stepper and propagator
  ConstantBField bField(Vector3D(0., 1_T, 0.));
  EigenStepper<
      ConstantBField,
      StepperExtensionList<DefaultExtension, DenseEnvironmentExtension>,
      detail::HighestValidAuctioneer>
      es(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DefaultExtension,
                                               DenseEnvironmentExtension>,
                          detail::HighestValidAuctioneer>,
             Navigator>
      prop(es, naviDet);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Manually set the extensions for each step and propagate through each
  // volume by propagation to the boundaries
  // Collect boundaries
  std::vector<Surface const*> surs;
  std::vector<std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>>
      boundaries = det->lowestTrackingVolume(tgContext, {0.5_m, 0., 0.})
                       ->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 1_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }
  boundaries =
      det->lowestTrackingVolume(tgContext, {1.5_m, 0., 0.})->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 2_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }
  boundaries =
      det->lowestTrackingVolume(tgContext, {2.5_m, 0., 0.})->boundarySurfaces();
  for (auto& b : boundaries) {
    if (b->surfaceRepresentation().center(tgContext).x() == 3_m) {
      surs.push_back(&(b->surfaceRepresentation()));
      break;
    }
  }

  // Build launcher through vacuum
  // Set options for propagator

  PropagatorOptions<ActionList<StepCollector, DebugOutput>,
                    AbortList<EndOfWorld>>
      propOptsDef(tgContext, mfContext);
  abortList.get<EndOfWorld>().maxX = 1_m;
  propOptsDef.abortList = abortList;
  propOptsDef.maxSteps = 100;
  propOptsDef.maxStepSize = 1.5_m;
  propOptsDef.debug = false;

  // Build stepper and propagator
  EigenStepper<ConstantBField, StepperExtensionList<DefaultExtension>> esDef(
      bField);
  Propagator<
      EigenStepper<ConstantBField, StepperExtensionList<DefaultExtension>>,
      Navigator>
      propDef(esDef, naviDet);

  // Launch and collect results
  const auto& resultDef =
      propDef.propagate(sbtp, *(surs[0]), propOptsDef).value();
  const StepCollector::this_result& stepResultDef =
      resultDef.get<typename StepCollector::result_type>();

  // Check the exit situation of the first volume
  std::pair<Vector3D, Vector3D> endParams, endParamsControl;
  for (unsigned int i = 0; i < stepResultDef.position.size(); i++) {
    if (1_m - stepResultDef.position[i].x() < 1e-4) {
      endParams =
          std::make_pair(stepResultDef.position[i], stepResultDef.momentum[i]);
      break;
    }
  }
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    if (1_m - stepResult.position[i].x() < 1e-4) {
      endParamsControl =
          std::make_pair(stepResult.position[i], stepResult.momentum[i]);
      break;
    }
  }

  if (propOptsDef.debug) {
    const auto debugString =
        resultDef.template get<DebugOutput::result_type>().debugString;
    std::cout << debugString << std::endl;
  }

  CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1_um);
  CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1_um);

  BOOST_TEST(endParams.first.x() == endParamsControl.first.x(),
             tt::tolerance(1e-5));
  BOOST_TEST(endParams.first.y() == endParamsControl.first.y(),
             tt::tolerance(1e-5));
  BOOST_TEST(endParams.first.z() == endParamsControl.first.z(),
             tt::tolerance(1e-5));
  BOOST_TEST(endParams.second.x() == endParamsControl.second.x(),
             tt::tolerance(1e-5));
  BOOST_TEST(endParams.second.y() == endParamsControl.second.y(),
             tt::tolerance(1e-5));
  BOOST_TEST(endParams.second.z() == endParamsControl.second.z(),
             tt::tolerance(1e-5));

  // Build launcher through material
  // Set initial parameters for the particle track by using the result of the
  // first volume
  startParams = endParams.first;
  startMom = endParams.second;
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtpPiecewise(
      cov, startParams, startMom, 1., 0.);

  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                AbortList<EndOfWorld>>
      propOptsDense(tgContext, mfContext);
  abortList.get<EndOfWorld>().maxX = 2_m;
  propOptsDense.abortList = abortList;
  propOptsDense.maxSteps = 1000;
  propOptsDense.maxStepSize = 1.5_m;
  propOptsDense.tolerance = 1e-8;

  // Build stepper and propagator
  EigenStepper<ConstantBField, StepperExtensionList<DenseEnvironmentExtension>>
      esDense(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DenseEnvironmentExtension>>,
             Navigator>
      propDense(esDense, naviDet);

  // Launch and collect results
  const auto& resultDense =
      propDense.propagate(sbtpPiecewise, *(surs[1]), propOptsDense).value();
  const StepCollector::this_result& stepResultDense =
      resultDense.get<typename StepCollector::result_type>();

  // Check the exit situation of the second volume
  for (unsigned int i = 0; i < stepResultDense.position.size(); i++) {
    if (2_m - stepResultDense.position[i].x() < 1e-4) {
      endParams = std::make_pair(stepResultDense.position[i],
                                 stepResultDense.momentum[i]);
      break;
    }
  }
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    if (2_m - stepResult.position[i].x() < 1e-4) {
      endParamsControl =
          std::make_pair(stepResult.position[i], stepResult.momentum[i]);
      break;
    }
  }

  CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1_um);
  CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1_um);
}

// Test case a). The DenseEnvironmentExtension should state that it is not
// valid in this case.
BOOST_AUTO_TEST_CASE(step_extension_trackercalomdt_test) {
  double rotationAngle = M_PI * 0.5;
  Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Vector3D yPos(0., 1., 0.);
  Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  MaterialProperties matProp(352.8, 407., 9.012, 4., 1.848e-3, 0.5_mm);

  CuboidVolumeBuilder cvb;
  CuboidVolumeBuilder::SurfaceConfig sConf1;
  sConf1.position = Vector3D(0.3_m, 0., 0.);
  sConf1.rotation.col(0) = xPos;
  sConf1.rotation.col(1) = yPos;
  sConf1.rotation.col(2) = zPos;
  sConf1.rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));
  sConf1.surMat = std::shared_ptr<const ISurfaceMaterial>(
      new HomogeneousSurfaceMaterial(matProp));
  sConf1.thickness = 1._mm;
  CuboidVolumeBuilder::LayerConfig lConf1;
  lConf1.surfaceCfg = sConf1;

  CuboidVolumeBuilder::SurfaceConfig sConf2;
  sConf2.position = Vector3D(0.6_m, 0., 0.);
  sConf2.rotation.col(0) = xPos;
  sConf2.rotation.col(1) = yPos;
  sConf2.rotation.col(2) = zPos;
  sConf2.rBounds =
      std::make_shared<const RectangleBounds>(RectangleBounds(0.5_m, 0.5_m));
  sConf2.surMat = std::shared_ptr<const ISurfaceMaterial>(
      new HomogeneousSurfaceMaterial(matProp));
  sConf2.thickness = 1._mm;
  CuboidVolumeBuilder::LayerConfig lConf2;
  lConf2.surfaceCfg = sConf2;

  CuboidVolumeBuilder::VolumeConfig muConf1;
  muConf1.position = {2.3_m, 0., 0.};
  muConf1.length = {20._cm, 20._cm, 20._cm};
  muConf1.volumeMaterial =
      std::shared_ptr<const IVolumeMaterial>(new HomogeneousVolumeMaterial(
          Material(352.8, 407., 9.012, 4., 1.848e-3)));
  muConf1.name = "MDT1";
  CuboidVolumeBuilder::VolumeConfig muConf2;
  muConf2.position = {2.7_m, 0., 0.};
  muConf2.length = {20._cm, 20._cm, 20._cm};
  muConf2.volumeMaterial =
      std::shared_ptr<const IVolumeMaterial>(new HomogeneousVolumeMaterial(
          Material(352.8, 407., 9.012, 4., 1.848e-3)));
  muConf2.name = "MDT2";

  CuboidVolumeBuilder::VolumeConfig vConf1;
  vConf1.position = {0.5_m, 0., 0.};
  vConf1.length = {1._m, 1._m, 1._m};
  vConf1.layerCfg = {lConf1, lConf2};
  vConf1.name = "Tracker";
  CuboidVolumeBuilder::VolumeConfig vConf2;
  vConf2.position = {1.5_m, 0., 0.};
  vConf2.length = {1._m, 1._m, 1._m};
  vConf2.volumeMaterial =
      std::shared_ptr<const IVolumeMaterial>(new HomogeneousVolumeMaterial(
          Material(352.8, 407., 9.012, 4., 1.848e-3)));
  vConf2.name = "Calorimeter";
  CuboidVolumeBuilder::VolumeConfig vConf3;
  vConf3.position = {2.5_m, 0., 0.};
  vConf3.length = {1._m, 1._m, 1._m};
  vConf3.volumeCfg = {muConf1, muConf2};
  vConf3.name = "Muon system";
  CuboidVolumeBuilder::Config conf;
  conf.volumeCfg = {vConf1, vConf2, vConf3};
  conf.position = {1.5_m, 0., 0.};
  conf.length = {3._m, 1._m, 1._m};

  // Build detector
  cvb.setConfig(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  // Build navigator
  Navigator naviVac(detector);
  naviVac.resolvePassive = true;
  naviVac.resolveMaterial = true;
  naviVac.resolveSensitive = true;

  // Set initial parameters for the particle track
  Covariance cov = Covariance::Identity();
  Vector3D startParams(0., 0., 0.), startMom(1._GeV, 0., 0.);
  SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(cov, startParams,
                                                       startMom, 1., 0.);

  // Set options for propagator
  DenseStepperPropagatorOptions<ActionList<StepCollector, MaterialInteractor>,
                                AbortList<EndOfWorld>>
      propOpts(tgContext, mfContext);
  propOpts.abortList.get<EndOfWorld>().maxX = 3._m;

  // Build stepper and propagator
  ConstantBField bField(Vector3D(0., 0., 0.));
  EigenStepper<
      ConstantBField,
      StepperExtensionList<DefaultExtension, DenseEnvironmentExtension>,
      detail::HighestValidAuctioneer>
      es(bField);
  Propagator<EigenStepper<ConstantBField,
                          StepperExtensionList<DefaultExtension,
                                               DenseEnvironmentExtension>,
                          detail::HighestValidAuctioneer>,
             Navigator>
      prop(es, naviVac);

  // Launch and collect results
  const auto& result = prop.propagate(sbtp, propOpts).value();
  const StepCollector::this_result& stepResult =
      result.get<typename StepCollector::result_type>();

  // Test that momentum changes only occured at the right detector parts
  double lastMomentum = stepResult.momentum[0].x();
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    // Test for changes
    if ((stepResult.position[i].x() > 0.3_m &&
         stepResult.position[i].x() < 0.6_m) ||
        (stepResult.position[i].x() > 0.6_m &&
         stepResult.position[i].x() <= 1._m) ||
        (stepResult.position[i].x() > 1._m &&
         stepResult.position[i].x() <= 2._m) ||
        (stepResult.position[i].x() > 2.2_m &&
         stepResult.position[i].x() <= 2.4_m) ||
        (stepResult.position[i].x() > 2.6_m &&
         stepResult.position[i].x() <= 2.8_m)) {
      BOOST_TEST(stepResult.momentum[i].x() <= lastMomentum);
      lastMomentum = stepResult.momentum[i].x();
    } else
    // Test the absence of momentum loss
    {
      if (stepResult.position[i].x() < 0.3_m ||
          (stepResult.position[i].x() > 2._m &&
           stepResult.position[i].x() <= 2.2_m) ||
          (stepResult.position[i].x() > 2.4_m &&
           stepResult.position[i].x() <= 2.6_m) ||
          (stepResult.position[i].x() > 2.8_m &&
           stepResult.position[i].x() <= 3._m)) {
        BOOST_TEST(stepResult.momentum[i].x() == lastMomentum);
      }
    }
  }
}
}  // namespace Test
}  // namespace Acts