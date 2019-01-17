// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Stepper Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include <fstream>

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"

// TODO: Testing of covariances in Integration test - requires N-layer box
// detector for implementation of DenseEnvironmentExtension
namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

  using cstep = detail::ConstrainedStep;

  ///
  /// @brief Aborter for the case that a particle leaves the detector or reaches
  /// a custom made threshold.
  ///
  struct EndOfWorld
  {
    /// Maximum value in x-direction of the detector
    double maxX = 1. * units::_m;

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
    bool
    operator()(propagator_state_t& state, const stepper_t& stepper) const
    {
      const double tolerance = state.options.targetTolerance;
      if (maxX - std::abs(stepper.position(state.stepping).x()) <= tolerance
          || std::abs(stepper.position(state.stepping).y()) >= 0.5 * units::_m
          || std::abs(stepper.position(state.stepping).z()) >= 0.5 * units::_m)
        return true;
      return false;
    }
  };

  ///
  /// @brief Data collector while propagation
  ///
  struct StepCollector
  {

    ///
    /// @brief Data container for result analysis
    ///
    struct this_result
    {
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
    void
    operator()(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      result.position.push_back(stepper.position(state.stepping));
      result.momentum.push_back(stepper.momentum(state.stepping)
                                * stepper.direction(state.stepping));
    }
  };

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
  BOOST_AUTO_TEST_CASE(step_extension_vacuum_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0.5 * units::_m, 0., 0.};
    vConf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {0.5 * units::_m, 0., 0.};
    conf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> vacuum = tgb.trackingGeometry();

    // Build navigator
    Navigator naviVac(vacuum);
    naviVac.resolvePassive   = true;
    naviVac.resolveMaterial  = true;
    naviVac.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList  = aList;
    propOpts.abortList   = abortList;
    propOpts.maxSteps    = 100;
    propOpts.maxStepSize = 0.5 * units::_m;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviVac);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Check that the propagation happend without interactions
    for (const auto& pos : stepResult.position) {
      CHECK_SMALL(pos.y(), 1. * units::_um);
      CHECK_SMALL(pos.z(), 1. * units::_um);
      if (pos == stepResult.position.back())
        CHECK_CLOSE_ABS(pos.x(), 1. * units::_m, 1. * units::_um);
    }
    for (const auto& mom : stepResult.momentum) {
      CHECK_CLOSE_ABS(mom, startMom, 1. * units::_keV);
    }

    // Rebuild and check the choice of extension
    ActionList<StepCollector> aListDef;

    // Set options for propagator
    PropagatorOptions<ActionList<StepCollector>, AbortList<EndOfWorld>>
        propOptsDef;
    propOptsDef.actionList  = aListDef;
    propOptsDef.abortList   = abortList;
    propOptsDef.maxSteps    = 100;
    propOptsDef.maxStepSize = 0.5 * units::_m;

    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension>>
        esDef(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension>>,
               Navigator>
        propDef(esDef, naviVac);

    // Launch and collect results
    const auto& resultDef = propDef.propagate(sbtp, propOptsDef);
    const StepCollector::this_result& stepResultDef
        = resultDef.get<typename StepCollector::result_type>();

    // Check that the right extension was chosen
    // If chosen correctly, the number of elements should be identical
    BOOST_TEST(stepResult.position.size() == stepResultDef.position.size());
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      CHECK_CLOSE_ABS(
          stepResult.position[i], stepResultDef.position[i], 1. * units::_um);
    }
    BOOST_TEST(stepResult.momentum.size() == stepResultDef.momentum.size());
    for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
      CHECK_CLOSE_ABS(
          stepResult.momentum[i], stepResultDef.momentum[i], 1. * units::_keV);
    }
  }
  // Test case b). The DefaultExtension should state that it is invalid here.
  BOOST_AUTO_TEST_CASE(step_extension_material_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0.5 * units::_m, 0., 0.};
    vConf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConf.material = std::make_shared<const Material>(
        Material(352.8, 394.133, 9.012, 4., 1.848e-3));
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {0.5 * units::_m, 0., 0.};
    conf.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> material = tgb.trackingGeometry();

    // Build navigator
    Navigator naviMat(material);
    naviMat.resolvePassive   = true;
    naviMat.resolveMaterial  = true;
    naviMat.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(5. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList  = aList;
    propOpts.abortList   = abortList;
    propOpts.maxSteps    = 100;
    propOpts.maxStepSize = 0.5 * units::_m;
    propOpts.debug       = true;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviMat);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Check that there occured interaction
    for (const auto& pos : stepResult.position) {
      CHECK_SMALL(pos.y(), 1. * units::_um);
      CHECK_SMALL(pos.z(), 1. * units::_um);
      if (pos == stepResult.position.front()) {
        CHECK_SMALL(pos.x(), 1. * units::_um);
      } else {
        BOOST_CHECK_GT(std::abs(pos.x()), 1. * units::_um);
      }
    }
    for (const auto& mom : stepResult.momentum) {
      CHECK_SMALL(mom.y(), 1. * units::_keV);
      CHECK_SMALL(mom.z(), 1. * units::_keV);
      if (mom == stepResult.momentum.front()) {
        CHECK_CLOSE_ABS(mom.x(), 5. * units::_GeV, 1. * units::_keV);
      } else {
        BOOST_CHECK_LT(mom.x(), 5. * units::_GeV);
      }
    }

    // Rebuild and check the choice of extension
    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOptsDense;
    propOptsDense.actionList  = aList;
    propOptsDense.abortList   = abortList;
    propOptsDense.maxSteps    = 100;
    propOptsDense.maxStepSize = 0.5 * units::_m;
    propOptsDense.debug       = true;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DenseEnvironmentExtension>>
        esDense(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DenseEnvironmentExtension>>,
               Navigator>
        propDense(esDense, naviMat);

    // Launch and collect results
    const auto& resultDense = propDense.propagate(sbtp, propOptsDense);
    const StepCollector::this_result& stepResultDense
        = resultDense.get<typename StepCollector::result_type>();

    // Check that the right extension was chosen
    // If chosen correctly, the number of elements should be identical
    BOOST_TEST(stepResult.position.size() == stepResultDense.position.size());
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      CHECK_CLOSE_ABS(
          stepResult.position[i], stepResultDense.position[i], 1. * units::_um);
    }
    BOOST_TEST(stepResult.momentum.size() == stepResultDense.momentum.size());
    for (unsigned int i = 0; i < stepResult.momentum.size(); i++) {
      CHECK_CLOSE_ABS(stepResult.momentum[i],
                      stepResultDense.momentum[i],
                      1. * units::_keV);
    }

    ////////////////////////////////////////////////////////////////////

    // Re-launch the configuration with magnetic field
    bField.setField(0., 1. * units::_T, 0.);
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        esB(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        propB(esB, naviMat);

    const auto& resultB = propB.propagate(sbtp, propOptsDense);
    const StepCollector::this_result& stepResultB
        = resultB.get<typename StepCollector::result_type>();

    // Check that there occured interaction
    for (const auto& pos : stepResultB.position) {
      if (pos == stepResultB.position.front()) {
        CHECK_SMALL(pos, 1. * units::_um);
      } else {
        BOOST_CHECK_GT(std::abs(pos.x()), 1. * units::_um);
        CHECK_SMALL(pos.y(), 1. * units::_um);
        BOOST_CHECK_GT(std::abs(pos.z()), 1. * units::_um);
      }
    }
    for (const auto& mom : stepResultB.momentum) {
      if (mom == stepResultB.momentum.front()) {
        CHECK_CLOSE_ABS(mom, startMom, 1. * units::_keV);
      } else {
        BOOST_CHECK_NE(mom.x(), 5. * units::_GeV);
        CHECK_SMALL(mom.y(), 1. * units::_keV);
        BOOST_CHECK_NE(mom.z(), 0.);
      }
    }
  }
  // Test case c). Both should be involved in their part of the detector
  BOOST_AUTO_TEST_CASE(step_extension_vacmatvac_test)
  {
    CuboidVolumeBuilder               cvb;
    CuboidVolumeBuilder::VolumeConfig vConfVac1;
    vConfVac1.position = {0.5 * units::_m, 0., 0.};
    vConfVac1.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConfVac1.name     = "First vacuum volume";
    CuboidVolumeBuilder::VolumeConfig vConfMat;
    vConfMat.position = {1.5 * units::_m, 0., 0.};
    vConfMat.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConfMat.material = std::make_shared<const Material>(
        Material(352.8, 394.133, 9.012, 4., 1.848e-3));
    vConfMat.name = "Material volume";
    CuboidVolumeBuilder::VolumeConfig vConfVac2;
    vConfVac2.position = {2.5 * units::_m, 0., 0.};
    vConfVac2.length   = {1. * units::_m, 1. * units::_m, 1. * units::_m};
    vConfVac2.name     = "Second vacuum volume";
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg = {vConfVac1, vConfMat, vConfVac2};
    conf.position  = {1.5 * units::_m, 0., 0.};
    conf.length    = {3. * units::_m, 1. * units::_m, 1. * units::_m};

    // Build detector
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        std::make_shared<const CuboidVolumeBuilder>(cvb));
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> det = tgb.trackingGeometry();

    // Build navigator
    Navigator naviDet(det);
    naviDet.resolvePassive   = true;
    naviDet.resolveMaterial  = true;
    naviDet.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov    = ActsSymMatrixD<5>::Identity();
    auto              covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(5. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector> aList;
    AbortList<EndOfWorld>     abortList;
    abortList.get<EndOfWorld>().maxX = 3. * units::_m;

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOpts;
    propOpts.actionList  = aList;
    propOpts.abortList   = abortList;
    propOpts.maxSteps    = 100;
    propOpts.maxStepSize = 0.5 * units::_m;

    // Build stepper and propagator
    ConstantBField bField(Vector3D(0., 1. * units::_T, 0.));
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension,
                                      DenseEnvironmentExtension>,
                 detail::HighestValidAuctioneer>
        es(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension,
                                                 DenseEnvironmentExtension>,
                            detail::HighestValidAuctioneer>,
               Navigator>
        prop(es, naviDet);

    // Launch and collect results
    const auto&                       result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    // Manually set the extensions for each step and propagate through each
    // volume by propagation to the boundaries
    // Collect boundaries
    std::vector<Surface const*> surs;
    std::vector<std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>>
        boundaries = det->lowestTrackingVolume({0.5 * units::_m, 0., 0.})
                         ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 1. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }
    boundaries = det->lowestTrackingVolume({1.5 * units::_m, 0., 0.})
                     ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 2. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }
    boundaries = det->lowestTrackingVolume({2.5 * units::_m, 0., 0.})
                     ->boundarySurfaces();
    for (auto& b : boundaries) {
      if (b->surfaceRepresentation().center().x() == 3. * units::_m) {
        surs.push_back(&(b->surfaceRepresentation()));
        break;
      }
    }

    // Build launcher through vacuum
    // Set options for propagator
    ActionList<StepCollector> aListDef;

    PropagatorOptions<ActionList<StepCollector>, AbortList<EndOfWorld>>
        propOptsDef;
    abortList.get<EndOfWorld>().maxX = 1. * units::_m;
    propOptsDef.actionList           = aListDef;
    propOptsDef.abortList            = abortList;
    propOptsDef.maxSteps             = 100;
    propOptsDef.maxStepSize          = 0.5 * units::_m;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DefaultExtension>>
        esDef(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DefaultExtension>>,
               Navigator>
        propDef(esDef, naviDet);

    // Launch and collect results
    const auto& resultDef = propDef.propagate(sbtp, *(surs[0]), propOptsDef);
    const StepCollector::this_result& stepResultDef
        = resultDef.get<typename StepCollector::result_type>();

    // Check the exit situation of the first volume
    std::pair<Vector3D, Vector3D> endParams, endParamsControl;
    for (unsigned int i = 0; i < stepResultDef.position.size(); i++) {
      if (1. * units::_m - stepResultDef.position[i].x() < 1e-4) {
        endParams = std::make_pair(stepResultDef.position[i],
                                   stepResultDef.momentum[i]);
        break;
      }
    }
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (1. * units::_m - stepResult.position[i].x() < 1e-4) {
        endParamsControl
            = std::make_pair(stepResult.position[i], stepResult.momentum[i]);
        break;
      }
    }

    CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1. * units::_um);
    CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1. * units::_um);

    // Build launcher through material
    // Set initial parameters for the particle track by using the result of the
    // first volume
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(
        ActsSymMatrixD<5>::Identity());
    startParams = endParams.first;
    startMom    = endParams.second;
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtpPiecewise(
        std::move(covPtr), startParams, startMom, 1.);

    // Set options for propagator
    DenseStepperPropagatorOptions<ActionList<StepCollector>,
                                  AbortList<EndOfWorld>>
        propOptsDense;
    abortList.get<EndOfWorld>().maxX = 2. * units::_m;
    propOptsDense.actionList         = aList;
    propOptsDense.abortList          = abortList;
    propOptsDense.maxSteps           = 100;
    propOptsDense.maxStepSize        = 0.5 * units::_m;
    propOptsDense.tolerance          = 1e-8;

    // Build stepper and propagator
    EigenStepper<ConstantBField,
                 VoidIntersectionCorrector,
                 StepperExtensionList<DenseEnvironmentExtension>>
        esDense(bField);
    Propagator<EigenStepper<ConstantBField,
                            VoidIntersectionCorrector,
                            StepperExtensionList<DenseEnvironmentExtension>>,
               Navigator>
        propDense(esDense, naviDet);

    // Launch and collect results
    const auto& resultDense
        = propDense.propagate(sbtpPiecewise, *(surs[1]), propOptsDense);
    const StepCollector::this_result& stepResultDense
        = resultDense.get<typename StepCollector::result_type>();

    // Check the exit situation of the second volume
    for (unsigned int i = 0; i < stepResultDense.position.size(); i++) {
      if (2. * units::_m - stepResultDense.position[i].x() < 1e-4) {
        endParams = std::make_pair(stepResultDense.position[i],
                                   stepResultDense.momentum[i]);
        break;
      }
    }
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      if (2. * units::_m - stepResult.position[i].x() < 1e-4) {
        endParamsControl
            = std::make_pair(stepResult.position[i], stepResult.momentum[i]);
        break;
      }
    }

    CHECK_CLOSE_ABS(endParams.first, endParamsControl.first, 1. * units::_um);
    CHECK_CLOSE_ABS(endParams.second, endParamsControl.second, 1. * units::_um);
  }
}  // namespace Test
}  // namespace Acts
