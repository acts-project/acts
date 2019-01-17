// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Propagator Tests
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

using Acts::VectorHelpers::perp;

namespace Acts {

namespace Test {

  using cstep = detail::ConstrainedStep;

  /// An observer that measures the perpendicular distance
  struct PerpendicularMeasure
  {

    /// Simple result struct to be returned
    struct this_result
    {
      double distance = std::numeric_limits<double>::max();
    };

    using result_type = this_result;

    PerpendicularMeasure() = default;

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      result.distance = perp(stepper.position(state.stepping));
    }

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& /*unused*/,
               const stepper_t& /*unused*/) const
    {
    }
  };

  /// An observer that measures the perpendicular distance
  template <typename Surface>
  struct SurfaceObserver
  {

    // the surface to be intersected
    const Surface* surface = nullptr;
    // the tolerance for intersection
    double tolerance = 1.e-5;

    /// Simple result struct to be returned
    struct this_result
    {
      size_t surfaces_passed  = 0;
      double surface_passed_r = std::numeric_limits<double>::max();
    };

    using result_type = this_result;

    SurfaceObserver() = default;

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      if (surface && !result.surfaces_passed) {
        // calculate the distance to the surface
        const double distance
            = surface
                  ->intersectionEstimate(stepper.position(state.stepping),
                                         stepper.direction(state.stepping),
                                         forward,
                                         true)
                  .pathLength;
        // Adjust the step size so that we cannot cross the target surface
        state.stepping.stepSize.update(distance, cstep::actor);
        // return true if you fall below tolerance
        if (std::abs(distance) <= tolerance) {
          ++result.surfaces_passed;
          result.surface_passed_r = perp(stepper.position(state.stepping));
          // release the step size, will be re-adjusted
          state.stepping.stepSize.release(cstep::actor);
        }
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& /*unused*/,
               const stepper_t& /*unused*/) const
    {
    }
  };

  /// An observer that scatters at a given pathlength
  /// with a fixed sigma(phi) and sigma(theta)
  struct PathScatterer
  {

    // the surface to be intersected
    double path_limit = std::numeric_limits<double>::max();
    // scattering deltas
    double sigma_phi   = 0.05;
    double sigma_theta = 0.05;
    // the tolerance for intersection
    double tolerance = 1.e-5;

    /// Simple result struct to be returned
    struct this_result
    {
      bool scattered = false;
    };

    using result_type = this_result;

    PathScatterer() = default;

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state,
               const stepper_t& /*unused*/,
               result_type& result) const
    {
      if (!result.scattered
          && std::abs(state.stepping.pathAccumulated - path_limit)
              < tolerance) {
        // now here we should apply the scattering
        result.scattered = true;
        // do the update and reinitialize the jacobians
        state.stepping.covarianceTransport(true);
        state.stepping.cov(ePHI, ePHI) += sigma_phi * sigma_phi;
        state.stepping.cov(eTHETA, eTHETA) += sigma_theta * sigma_theta;
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& /*unused*/,
               const stepper_t& /*unused*/) const
    {
    }
  };

  // Global definitions
  // The path limit abort
  using path_limit = detail::PathLimitReached;

  using BField_type         = ConstantBField;
  using EigenStepper_type   = EigenStepper<BField_type>;
  using EigenPropagatorType = Propagator<EigenStepper_type>;

  const double        Bz = 2. * units::_T;
  BField_type         bField(0, 0, Bz);
  EigenStepper_type   estepper(bField);
  EigenPropagatorType epropagator(std::move(estepper));

  auto mSurface
      = Surface::makeShared<CylinderSurface>(nullptr, 10., 1000. * units::_mm);
  auto cSurface
      = Surface::makeShared<CylinderSurface>(nullptr, 150., 1000. * units::_mm);

  const int ntests = 5;

  // This tests the Options
  BOOST_AUTO_TEST_CASE(PropagatorOptions_)
  {

    using null_options_type = PropagatorOptions<>;
    null_options_type null_options;
    // todo write null options test

    using ActionList_type      = ActionList<PerpendicularMeasure>;
    using AbortConditions_type = AbortList<>;

    using options_type
        = PropagatorOptions<ActionList_type, AbortConditions_type>;

    options_type options;
  }

  BOOST_DATA_TEST_CASE(
      cylinder_passage_observer_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    using CylinderObserver     = SurfaceObserver<CylinderSurface>;
    using ActionList_type      = ActionList<CylinderObserver>;
    using AbortConditions_type = AbortList<>;

    // setup propagation options
    PropagatorOptions<ActionList_type, AbortConditions_type> options;

    options.pathLimit   = 20 * units::_m;
    options.maxStepSize = 1 * units::_cm;

    // set the surface to be passed by
    options.actionList.get<CylinderObserver>().surface = mSurface.get();

    using so_result = typename CylinderObserver::result_type;

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = dcharge;
    Vector3D              pos(x, y, z);
    Vector3D              mom(px, py, pz);
    CurvilinearParameters start(nullptr, pos, mom, q);
    // propagate to the cylinder surface
    const auto& result = epropagator.propagate(start, *cSurface, options);
    auto&       sor    = result.get<so_result>();

    BOOST_CHECK_EQUAL(sor.surfaces_passed, 1);
    CHECK_CLOSE_ABS(sor.surface_passed_r, 10., 1e-5);
  }

  BOOST_DATA_TEST_CASE(
      curvilinear_additive_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    // setup propagation options - the tow step options
    PropagatorOptions<> options_2s;
    options_2s.pathLimit   = 50 * units::_cm;
    options_2s.maxStepSize = 1 * units::_cm;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);
    // propagate to a path length of 100 with two steps of 50
    const auto& mid_parameters
        = epropagator.propagate(start, options_2s).endParameters;
    const auto& end_parameters_2s
        = epropagator.propagate(*mid_parameters, options_2s).endParameters;

    // setup propagation options - the one step options
    PropagatorOptions<> options_1s;
    options_1s.pathLimit   = 100 * units::_cm;
    options_1s.maxStepSize = 1 * units::_cm;
    // propagate to a path length of 100 in one step
    const auto& end_parameters_1s
        = epropagator.propagate(start, options_1s).endParameters;

    // test that the propagation is additive
    CHECK_CLOSE_REL(
        end_parameters_1s->position(), end_parameters_2s->position(), 0.001);

    const auto& cov_1s = *(end_parameters_1s->covariance());
    const auto& cov_2s = *(end_parameters_2s->covariance());

    CHECK_CLOSE_COVARIANCE(cov_1s, cov_2s, 0.001);
  }

  BOOST_DATA_TEST_CASE(
      cylinder_additive_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    // setup propagation options - 2 setp options
    PropagatorOptions<> options_2s;
    options_2s.pathLimit   = 10 * units::_m;
    options_2s.maxStepSize = 1 * units::_cm;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);
    // propagate to a final surface with one stop in between
    const auto& mid_parameters
        = epropagator.propagate(start, *mSurface, options_2s).endParameters;

    const auto& end_parameters_2s
        = epropagator.propagate(*mid_parameters, *cSurface, options_2s)
              .endParameters;

    // setup propagation options - one step options
    PropagatorOptions<> options_1s;
    options_1s.pathLimit   = 10 * units::_m;
    options_1s.maxStepSize = 1 * units::_cm;
    // propagate to a final surface in one stop
    const auto& end_parameters_1s
        = epropagator.propagate(start, *cSurface, options_1s).endParameters;

    // test that the propagation is additive
    CHECK_CLOSE_REL(
        end_parameters_1s->position(), end_parameters_2s->position(), 0.001);

    const auto& cov_1s = *(end_parameters_1s->covariance());
    const auto& cov_2s = *(end_parameters_2s->covariance());

    CHECK_CLOSE_COVARIANCE(cov_1s, cov_2s, 0.001);
  }

}  // namespace Test
}  // namespace Acts
