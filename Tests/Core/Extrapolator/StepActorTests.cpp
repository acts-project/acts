// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE StepActor Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/StepActor.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "DetectorBuilder.hpp"

namespace Acts {
namespace Test {

  ///
  /// @brief Aborter for the case that a particle leaves the detector
  ///
  struct EndOfWorld
  {
    /// @brief Constructor
    EndOfWorld() = default;

    /// @brief Main call operator for the abort operation
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    /// @return Boolean statement if the particle is still in the detector
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (std::abs(state.stepping.position().x()) >= 2. * units::_m
          || std::abs(state.stepping.position().y()) >= 0.5 * units::_m
          || std::abs(state.stepping.position().z()) >= 0.5 * units::_m)
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
      // Covariance matrix of the propagator after each step
      std::vector<ActsSymMatrixD<5>> cov;
    };

    using result_type = this_result;

    /// @brief Main call operator for the action list. It stores the data for
    /// analysis afterwards
    ///
    /// @tparam propagator_state_t Type of the propagator state
    /// @param [in] state State of the propagator
    /// @param [out] result Struct which is filled with the data
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      result.position.push_back(state.stepping.position());
      result.momentum.push_back(state.stepping.momentum());
      result.cov.push_back(state.stepping.cov);
    }
  };

  /// @brief This is a test for the StepActor that updates the propagator in
  /// dense material. The test consists out of 2 parts.
  /// 1.: Test that the StepActor does not act if there is no material in a
  /// TrackingVolume
  /// 2.: Test that the StepActor acts if there is material in a TrackingVolume
  BOOST_AUTO_TEST_CASE(step_actor_test)
  {
    {
      // Build detector
      std::shared_ptr<TrackingGeometry> vacuum = buildVacDetector();

      // Build navigator
      Navigator naviVac(vacuum);
      naviVac.resolvePassive   = true;
      naviVac.resolveMaterial  = true;
      naviVac.resolveSensitive = true;

      // Set initial parameters for the particle track
      ActsSymMatrixD<5> cov;
      cov << 1. * units::_mm, 0., 0., 0., 0., 0., 1. * units::_mm, 0., 0., 0.,
          0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.;
      auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
      Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
      SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
          std::move(covPtr), startParams, startMom, 1.);

      // Create action list for surface collection
      ActionList<StepCollector, StepActor> aList;
      AbortList<EndOfWorld> abortList;

      // Set options for propagator
      Propagator<EigenStepper<ConstantBField>, Navigator>::
          Options<ActionList<StepCollector, StepActor>, AbortList<EndOfWorld>>
              propOpts;
      propOpts.actionList     = aList;
      propOpts.stopConditions = abortList;
      propOpts.maxSteps       = 1e6;
      propOpts.maxStepSize    = 2. * units::_m;

      // Re-configure propagation with B-field
      ConstantBField               bField(Vector3D(0., 0., 0.));
      EigenStepper<ConstantBField> es(bField);
      Propagator<EigenStepper<ConstantBField>, Navigator> prop(es, naviVac);

      // Launch and collect results
      const auto&                       result = prop.propagate(sbtp, propOpts);
      const StepCollector::this_result& stepResult
          = result.get<typename StepCollector::result_type>();

      // Check that the propagation happend in a straight line without
      // interactions
      for (const auto& pos : stepResult.position) {
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() == 0.);
        if (pos == stepResult.position.back())
          BOOST_TEST(pos.x() == 2. * units::_m);
      }
      for (const auto& mom : stepResult.momentum) {
        BOOST_TEST(mom.x() == 1. * units::_GeV);
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() == 0.);
      }
      for (const auto& c : stepResult.cov) {
        BOOST_TEST(c == ActsSymMatrixD<5>::Identity());
      }
    }
    {
      // Build detector
      std::shared_ptr<TrackingGeometry> material = buildMatDetector();

      // Build navigator
      Navigator naviMat(material);
      naviMat.resolvePassive   = true;
      naviMat.resolveMaterial  = true;
      naviMat.resolveSensitive = true;

      // Set initial parameters for the particle track
      ActsSymMatrixD<5> cov;
      cov << 1. * units::_mm, 0., 0., 0., 0., 0., 1. * units::_mm, 0., 0., 0.,
          0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.;
      auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
      Vector3D startParams(0., 0., 0.), startMom(1. * units::_GeV, 0., 0.);
      SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
          std::move(covPtr), startParams, startMom, 1.);

      // Create action list for surface collection
      ActionList<StepCollector, StepActor> aList;
      AbortList<EndOfWorld> abortList;

      // Set options for propagator
      Propagator<EigenStepper<ConstantBField>, Navigator>::
          Options<ActionList<StepCollector, StepActor>, AbortList<EndOfWorld>>
              propOpts;
      propOpts.actionList     = aList;
      propOpts.stopConditions = abortList;
      propOpts.maxSteps       = 1e6;
      propOpts.maxStepSize    = 2. * units::_m;
      propOpts.debug          = true;

      // Re-configure propagation with B-field
      ConstantBField               bField(Vector3D(0., 0., 0.));
      EigenStepper<ConstantBField> es(bField);
      Propagator<EigenStepper<ConstantBField>, Navigator> prop(es, naviMat);

      // Launch and collect results
      const auto&                       result = prop.propagate(sbtp, propOpts);
      const StepCollector::this_result& stepResult
          = result.get<typename StepCollector::result_type>();

      // Check that there occured interaction
      for (const auto& pos : stepResult.position) {
        BOOST_TEST(pos.y() == 0.);
        BOOST_TEST(pos.z() == 0.);
        if (pos == stepResult.position.front()) {
          BOOST_TEST(pos.x() == 0.);
        } else {
          BOOST_TEST(pos.x() != 0.);
        }
      }
      for (const auto& mom : stepResult.momentum) {
        BOOST_TEST(mom.y() == 0.);
        BOOST_TEST(mom.z() == 0.);
        if (mom == stepResult.momentum.front()) {
          BOOST_TEST(mom.x() == 1. * units::_GeV);
        } else {
          BOOST_TEST(mom.x() != 1. * units::_GeV);
        }
      }
      for (const auto& c : stepResult.cov) {
        if (c == stepResult.cov.front()) {
          BOOST_TEST(c == ActsSymMatrixD<5>::Identity());
        } else {
          BOOST_TEST(c != ActsSymMatrixD<5>::Identity());
        }
      }
      //////////////////////////////////////////////////////////////////

      // Re-launch the configuration with magnetic field
      bField.setField(0., 1. * units::_T, 0.);
      EigenStepper<ConstantBField> esB(bField);
      Propagator<EigenStepper<ConstantBField>, Navigator> probB(esB, naviMat);

      const auto& resultB = probB.propagate(sbtp, propOpts);
      const StepCollector::this_result& stepResultB
          = resultB.get<typename StepCollector::result_type>();

      // Check that there occured interaction
      for (const auto& pos : stepResultB.position) {
        if (pos == stepResultB.position.front()) {
          BOOST_TEST(pos.x() == 0.);
          BOOST_TEST(pos.y() == 0.);
          BOOST_TEST(pos.z() == 0.);
        } else {
          BOOST_TEST(pos.x() != 0.);
          BOOST_TEST(pos.y() == 0.);
          BOOST_TEST(pos.z() != 0.);
        }
      }
      for (const auto& mom : stepResultB.momentum) {
        if (mom == stepResultB.momentum.front()) {
          BOOST_TEST(mom.x() == 1. * units::_GeV);
          BOOST_TEST(mom.y() == 0.);
          BOOST_TEST(mom.z() == 0.);
        } else {
          BOOST_TEST(mom.x() != 1. * units::_GeV);
          BOOST_TEST(mom.y() == 0.);
          BOOST_TEST(mom.z() != 0.);
        }
      }
      for (const auto& c : stepResultB.cov) {
        if (c == stepResultB.cov.front()) {
          BOOST_TEST(c == ActsSymMatrixD<5>::Identity());
        } else {
          BOOST_TEST(c != ActsSymMatrixD<5>::Identity());
        }
      }
    }
  }
  // TODO: Test step size release

}  // namespace Test
}  // namespace Acts
