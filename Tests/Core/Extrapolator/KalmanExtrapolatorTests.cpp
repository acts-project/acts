// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Kalman Extrapolator Test
#include <boost/test/included/unit_test.hpp>

#include <math.h>
#include <random>
#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/Common/CubicTrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

  using Jacobian = ActsMatrixD<5, 5>;

  ///
  /// @brief the bound state propagation
  ///
  struct StepWiseActor
  {

    /// The result is the piece-wise jacobian
    struct this_result
    {

      std::vector<Jacobian> jacobians = {};
      std::vector<double>   paths     = {};

      Jacobian fullJacobian = Jacobian::Identity();
      double   fullPath     = 0.;

      bool finalized = false;
    };
    /// Broadcast the result type
    using result_type = this_result;

    /// @brief Kalman sequence operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    ///
    /// @param state is the mutable propagator state object
    /// @param result is the mutable result state object
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& state, result_type& result) const
    {
      // Listen to the surface and create bound state where necessary
      auto surface = state.navigation.currentSurface;
      if (surface and surface->associatedDetectorElement()) {
        // Create a bound state and log the jacobian
        auto boundState = state.stepping.boundState(*surface, true);
        result.jacobians.push_back(std::move(std::get<Jacobian>(boundState)));
        result.paths.push_back(std::get<double>(boundState));
      }
      // Also store the
      if ((state.navigation.navigationBreak or state.navigation.targetReached)
          and not result.finalized) {
        // Set the last stepping parameter
        result.paths.push_back(state.stepping.pathAccumulated);
        // Set the full parameters
        result.fullJacobian = state.stepping.jacobian;
        result.fullPath     = state.stepping.pathAccumulated;
        // Remember that you finalized this
        result.finalized = true;
      }
    }

    /// @brief Kalman sequence operation - void operation
    ///
    /// @tparam propagator_state_t is the type of Propagagor state
    ///
    /// @param state is the mutable propagator state object
    template <typename propagator_state_t>
    void
    operator()(propagator_state_t& /*state*/) const
    {
    }
  };

  ///
  /// @brief Unit test for Kalman fitter propagation
  ///
  BOOST_AUTO_TEST_CASE(kalman_extrapolator)
  {
    // Build detector, take the cubic detector
    CubicTrackingGeometry cGeometry;
    auto                  detector = cGeometry();

    // The Naivgator through the detector geometry
    Navigator navigator(detector);
    navigator.resolvePassive   = false;
    navigator.resolveMaterial  = true;
    navigator.resolveSensitive = true;

    // Configure propagation with deactivated B-field
    ConstantBField bField(Vector3D(0., 0., 0.));
    using Stepper = EigenStepper<ConstantBField>;
    Stepper stepper(bField);
    using Propagator = Propagator<Stepper, Navigator>;
    Propagator propagator(stepper, navigator);

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    // The start position and start parameters
    Vector3D pos(-3. * units::_m, 0., 0.), mom(1. * units::_GeV, 0., 0);
    SingleCurvilinearTrackParameters<ChargedPolicy> start(
        std::move(covPtr), pos, mom, 1.);

    // Create the ActionList and AbortList
    using StepWiseResult = StepWiseActor::result_type;
    using StepWiseActors = ActionList<StepWiseActor>;
    using Aborters       = AbortList<detail::EndOfWorldReached>;

    // Create some options
    using StepWiseOptions =
        PropagatorOptions<StepWiseActors, Aborters>;
    StepWiseOptions swOptions;

    using PlainActors = ActionList<>;
    using PlainOptions =
        PropagatorOptions<PlainActors, Aborters>;
    PlainOptions pOptions;

    // Run the standard propagation
    const auto& pResult = propagator.propagate(start, pOptions);
    // Let's get the end parameters and jacobian matrix
    const auto& pJacobian = *(pResult.transportJacobian);

    // Run the stepwise propagation
    const auto& swResult       = propagator.propagate(start, swOptions);
    auto        swJacobianTest = swResult.template get<StepWiseResult>();

    // (1) Path length test
    double accPath = 0.;
    auto   swPaths = swJacobianTest.paths;
    // Sum up the step-wise paths, they are total though
    for (auto cpath = swPaths.rbegin(); cpath != swPaths.rend(); ++cpath) {
      if (cpath != swPaths.rend() - 1) {
        accPath += (*cpath) - (*(cpath + 1));
        continue;
      }
      accPath += (*cpath) - 0.;
    }
    BOOST_CHECK_EQUAL(swJacobianTest.fullPath, accPath);

    // (2) Jacobian test
    Jacobian accJacobian = Jacobian::Identity();
    // The stepwise jacobians
    auto swJacobians = swJacobianTest.jacobians;
    // The last-step jacobian, needed for the full jacobian transport
    const auto& swlJacobian = *(swResult.transportJacobian);

    // Build up the step-wise jacobians
    for (auto& j : swJacobians) {
      accJacobian = j * accJacobian;
    }
    accJacobian = swlJacobian * accJacobian;
    BOOST_TEST(pJacobian.isApprox(accJacobian));
  }
}
}