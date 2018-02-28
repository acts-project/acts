// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Propagator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/ActionList.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Propagator/detail/standard_abort_conditions.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using namespace propagation;

namespace Test {

  /// An observer that measures the perpendicular distance
  struct PerpendicularMeasure
  {

    /// Simple result struct to be returned
    struct this_result
    {
      double distance = std::numeric_limits<double>::max();
    };

    typedef this_result result_type;

    PerpendicularMeasure() {}

    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      result.distance = cache.pos.perp();
    }

    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
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

    typedef this_result result_type;

    SurfaceObserver() {}

    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      if (surface && !result.surfaces_passed) {
        // calculate the distance to the surface
        const double distance
            = surface
                  ->intersectionEstimate(
                      cache.position(), cache.direction(), true, true)
                  .pathLength;
        // Adjust the step size so that we cannot cross the target surface
        if (std::abs(cache.step_size) > std::abs(distance))
          cache.step_size = distance;
        // return true if you fall below tolerance
        if (std::abs(distance) <= tolerance) {
          ++result.surfaces_passed;
          result.surface_passed_r = cache.position().perp();
        }
      }
    }

    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
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

    typedef this_result result_type;

    PathScatterer() {}

    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      if (!result.scattered
          && std::abs(cache.accumulated_path - path_limit) < tolerance) {
        // now here we should apply the scattering
        result.scattered = true;
        // do the update and reinitialize the jacobians
        cache.apply_cov_transport(true);
        cache.cov(ePHI, ePHI) += sigma_phi * sigma_phi;
        cache.cov(eTHETA, eTHETA) += sigma_theta * sigma_theta;
      }
    }

    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
    }
  };

  /// An observer that scatters at a given surface
  /// with a fixed sigma(phi) and sigma(theta)
  template <typename Surface>
  struct SurfaceScatterer
  {

    // the surface to be intersected
    const Surface* surface = nullptr;
    // scattering deltas
    double sigma_phi   = 0.05;
    double sigma_theta = 0.05;
    // the tolerance for intersection
    double tolerance = 1.e-5;

    /// Simple result struct to be returned
    struct this_result
    {
      bool surfaces_passed = false;
    };

    typedef this_result result_type;

    SurfaceScatterer() {}

    template <typename cache_t>
    void
    operator()(cache_t& cache, result_type& result) const
    {
      if (surface && !result.surfaces_passed) {
        // calculate the distance to the surface
        const double distance
            = surface
                  ->intersectionEstimate(
                      cache.position(), cache.direction(), true, true)
                  .pathLength;
        // Adjust the step size so that we cannot cross the target surface
        if (std::abs(cache.step_size) > std::abs(distance))
          cache.step_size = distance;
        // return true if you fall below tolerance
        if (std::abs(distance) <= tolerance) {
          result.surfaces_passed = true;
          // do the update and reinitialize the jacobians
          cache.apply_cov_transport(*surface, true);
          cache.cov(ePHI, ePHI) += sigma_phi * sigma_phi;
          cache.cov(eTHETA, eTHETA) += sigma_theta * sigma_theta;
          // @HACK - release the step size - this should be triggered
          cache.step_size = 10. * units::_mm;
        }
      }
    }

    template <typename cache_t>
    void
    operator()(cache_t& cache) const
    {
      (void)cache;
    }
  };

  // Global definitions
  // The path limit abort
  typedef detail::path_limit_reached path_limit;

  typedef ConstantBField                BField_type;
  typedef EigenStepper<BField_type>     EigenStepper_type;
  typedef Propagator<EigenStepper_type> EigenPropagator_type;

  const double         Bz = 2. * units::_T;
  BField_type          bField(0, 0, Bz);
  EigenStepper_type    estepper(bField);
  EigenPropagator_type epropagator(std::move(estepper));

  CylinderSurface mSurface(nullptr, 10., 1000. * units::_mm);
  CylinderSurface cSurface(nullptr, 150., 1000. * units::_mm);

  const int ntests = 100;

  // This tests the Options
  BOOST_AUTO_TEST_CASE(PropgatorOptions_)
  {

    typedef typename Propagator<EigenStepper_type>::template Options<>
                      null_options_type;
    null_options_type null_options;
    // todo write null options test

    typedef ActionList<PerpendicularMeasure> ObsList_type;
    typedef AbortList<>                      AbortConditions_type;

    typedef typename Propagator<EigenStepper_type>::
        template Options<ObsList_type, AbortConditions_type>
            options_type;

    options_type options;
  }

  BOOST_DATA_TEST_CASE(
      cylinder_intersection_observer_,
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

    typedef ActionList<PerpendicularMeasure> ObsList_type;
    typedef AbortList<>                      AbortConditions_type;

    // setup propagation options
    typename EigenPropagator_type::template Options<ObsList_type,
                                                    AbortConditions_type>
        options;

    options.max_path_length = 20 * units::_m;
    options.max_step_size   = 1 * units::_cm;

    typedef typename PerpendicularMeasure::result_type pm_result;

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
    const auto& result = epropagator.propagate(start, cSurface, options);
    auto&       pmr    = result.get<pm_result>();

    // Test the end position
    BOOST_TEST(std::abs(result.endParameters->position().perp() - 150.)
               < 10e-4);
    BOOST_TEST(std::abs(pmr.distance - 150.) < 10e-4);
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

    typedef SurfaceObserver<CylinderSurface> CylinderObserver;
    typedef ActionList<CylinderObserver>     ObsList_type;
    typedef AbortList<>                      AbortConditions_type;

    // setup propagation options
    typename EigenPropagator_type::template Options<ObsList_type,
                                                    AbortConditions_type>
        options;

    options.max_path_length = 20 * units::_m;
    options.max_step_size   = 1 * units::_cm;

    // set the surface to be passed by
    options.action_list.get<CylinderObserver>().surface = &mSurface;

    typedef typename CylinderObserver::result_type so_result;

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
    const auto& result = epropagator.propagate(start, cSurface, options);
    auto&       sor    = result.get<so_result>();

    BOOST_TEST(sor.surfaces_passed == 1);
    BOOST_TEST(std::abs(sor.surface_passed_r - 10.) < 1e-5);
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
    typename EigenPropagator_type::template Options<> options_2s;
    options_2s.max_path_length = 50 * units::_cm;
    options_2s.max_step_size   = 1 * units::_cm;

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
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);
    // propagate to a path length of 100 with two steps of 50
    const auto& mid_parameters
        = epropagator.propagate(start, options_2s).endParameters;
    const auto& end_parameters_2s
        = epropagator.propagate(*mid_parameters, options_2s).endParameters;

    // setup propagation options - the one step options
    typename EigenPropagator_type::template Options<> options_1s;
    options_1s.max_path_length = 100 * units::_cm;
    options_1s.max_step_size   = 1 * units::_cm;
    // propagate to a path length of 100 in one step
    const auto& end_parameters_1s
        = epropagator.propagate(start, options_1s).endParameters;

    // test that the propagation is additive
    BOOST_TEST(end_parameters_1s->position().isApprox(
        end_parameters_2s->position(), 0.001));

    const auto& cov_1s = *(end_parameters_1s->covariance());
    const auto& cov_2s = *(end_parameters_2s->covariance());

    BOOST_TEST(cov_1s.isApprox(cov_2s, 0.001));
  }

  BOOST_DATA_TEST_CASE(
      curvilinear_scattering_,
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
          ^ bdata::random((bdata::seed = 4,
                           bdata::distribution
                           = std::normal_distribution<>(0, 0.05)))
          ^ bdata::random((bdata::seed = 5,
                           bdata::distribution
                           = std::normal_distribution<>(0, 0.05)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      sigma_phi,
      sigma_theta,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    // setup propagation options - the tow step options
    typename EigenPropagator_type::template Options<> options_2s;
    options_2s.max_path_length = 50 * units::_cm;
    options_2s.max_step_size   = 1 * units::_cm;

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
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);
    // propagate to a path length of 100 with two steps of 50
    const auto& mid_parameters
        = epropagator.propagate(start, options_2s).endParameters;
    // now 'scatter' in between and
    ActsSymMatrixD<5> mid_cov(*(mid_parameters->covariance()));
    mid_cov(ePHI, ePHI) += sigma_phi * sigma_phi;
    mid_cov(eTHETA, eTHETA) += sigma_theta * sigma_theta;

    auto mid_pos     = mid_parameters->position();
    auto mid_mom     = mid_parameters->momentum();
    auto mid_cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(mid_cov);
    CurvilinearParameters mid_pararameters_s(
        std::move(mid_cov_ptr), mid_pos, mid_mom, q);
    // the two step parameters now include the scattering ad the mid point
    const auto& end_parameters_2s_s
        = epropagator.propagate(mid_pararameters_s, options_2s).endParameters;

    // setup propagation options - the one step options
    typename EigenPropagator_type::template Options<> options_1s;
    options_1s.max_path_length = 100 * units::_cm;
    options_1s.max_step_size   = 1 * units::_cm;
    // propagate to a path length of 100 in one step
    const auto& end_parameters_1s
        = epropagator.propagate(start, options_1s).endParameters;

    // test that the propagation is additive
    // positions should still be the same
    BOOST_TEST(end_parameters_1s->position().isApprox(
        end_parameters_2s_s->position(), 0.001));

    const auto& cov_1s   = *(end_parameters_1s->covariance());
    const auto& cov_2s_s = *(end_parameters_2s_s->covariance());

    // covariances need to be different - check for difference here
    BOOST_TEST(!cov_1s.isApprox(cov_2s_s));

    // to get the covariances agree, we need an Observer that does the
    // scattering
    typedef ActionList<PathScatterer> ObsList_type;
    typedef AbortList<>               AbortConditions_type;

    // setup propagation options - the 1 step with scattering optios
    typename EigenPropagator_type::template Options<ObsList_type,
                                                    AbortConditions_type>
          options_1s_s;
    auto& _s       = options_1s_s.action_list.get<PathScatterer>();
    _s.path_limit  = 50. * units::_cm;
    _s.sigma_phi   = sigma_phi;
    _s.sigma_theta = sigma_theta;

    options_1s_s.max_path_length = 100 * units::_cm;
    options_1s_s.max_step_size   = 1 * units::_cm;
    // propagate to a path length of 100 in one step
    const auto& end_parameters_1s_s
        = epropagator.propagate(start, options_1s_s).endParameters;

    const auto& cov_1s_s = *(end_parameters_1s_s->covariance());

    // now it should be the same again, the PathScatter did the trick
    BOOST_TEST(cov_1s_s.isApprox(cov_2s_s));
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
    typename EigenPropagator_type::template Options<> options_2s;
    options_2s.max_path_length = 10 * units::_m;
    options_2s.max_step_size   = 1 * units::_cm;

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
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);
    // propagate to a final surface with one stop in between
    const auto& mid_parameters
        = epropagator.propagate(start, mSurface, options_2s).endParameters;

    const auto& end_parameters_2s
        = epropagator.propagate(*mid_parameters, cSurface, options_2s)
              .endParameters;

    // setup propagation options - one step options
    typename EigenPropagator_type::template Options<> options_1s;
    options_1s.max_path_length = 10 * units::_m;
    options_1s.max_step_size   = 1 * units::_cm;
    // propagate to a final surface in one stop
    const auto& end_parameters_1s
        = epropagator.propagate(start, cSurface, options_1s).endParameters;

    // test that the propagation is additive
    BOOST_TEST(end_parameters_1s->position().isApprox(
        end_parameters_2s->position(), 0.001));

    const auto& cov_1s = *(end_parameters_1s->covariance());
    const auto& cov_2s = *(end_parameters_2s->covariance());

    BOOST_TEST(cov_1s.isApprox(cov_2s, 0.001));
  }

  BOOST_DATA_TEST_CASE(
      cylinder_scattering_,
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
          ^ bdata::random((bdata::seed = 4,
                           bdata::distribution
                           = std::normal_distribution<>(0, 0.05)))
          ^ bdata::random((bdata::seed = 5,
                           bdata::distribution
                           = std::normal_distribution<>(0, 0.05)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      sigma_phi,
      sigma_theta,
      index)
  {
    double dcharge = -1 + 2 * charge;
    (void)index;

    // setup propagation options - 2 setp options
    typename EigenPropagator_type::template Options<> options_2s;
    options_2s.max_path_length = 10 * units::_m;
    options_2s.max_step_size   = 1 * units::_cm;

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
    auto cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(cov_ptr), pos, mom, q);

    // propagate to a final surface with one stop in between
    const auto& mid_parameters
        = epropagator.propagate(start, mSurface, options_2s).endParameters;
    // now 'scatter' on the mid surface
    ActsSymMatrixD<5> mid_cov(*(mid_parameters->covariance()));
    mid_cov(ePHI, ePHI) += sigma_phi * sigma_phi;
    mid_cov(eTHETA, eTHETA) += sigma_theta * sigma_theta;

    auto mid_pos     = mid_parameters->position();
    auto mid_mom     = mid_parameters->momentum();
    auto mid_cov_ptr = std::make_unique<const ActsSymMatrixD<5>>(mid_cov);
    BoundParameters mid_pararameters_s(
        std::move(mid_cov_ptr), mid_pos, mid_mom, q, mSurface);
    // the end parameters now carry the scattring of the mid surface
    const auto& end_parameters_2s_s
        = epropagator.propagate(mid_pararameters_s, cSurface, options_2s)
              .endParameters;

    // setup propagation options - one step options
    typename EigenPropagator_type::template Options<> options_1s;
    options_1s.max_path_length = 10 * units::_m;
    options_1s.max_step_size   = 1 * units::_cm;
    // propagate to a final surface in one stop
    const auto& end_parameters_1s
        = epropagator.propagate(start, cSurface, options_1s).endParameters;

    // test that the propagation is additive
    BOOST_TEST(end_parameters_1s->position().isApprox(
        end_parameters_2s_s->position(), 0.001));

    const auto& cov_1s   = *(end_parameters_1s->covariance());
    const auto& cov_2s_s = *(end_parameters_2s_s->covariance());

    BOOST_TEST(!cov_1s.isApprox(cov_2s_s, 0.001));

    // Cylinder Surface Scatterer
    typedef SurfaceScatterer<CylinderSurface> CylinderScatterer;

    // to get the covariances agree, we need an Observer that does the
    // scattering
    typedef ActionList<CylinderScatterer> ObsList_type;
    typedef AbortList<>                   AbortConditions_type;

    // setup propagation options - the 1 step with scattering optios
    typename EigenPropagator_type::template Options<ObsList_type,
                                                    AbortConditions_type>
        options_1s_s;
    options_1s_s.max_step_size = 1 * units::_cm;

    auto& _s       = options_1s_s.action_list.get<CylinderScatterer>();
    _s.surface     = &mSurface;
    _s.sigma_phi   = sigma_phi;
    _s.sigma_theta = sigma_theta;
    // propagate to a path length of 100 in one step
    const auto& end_parameters_1s_s
        = epropagator.propagate(start, cSurface, options_1s_s).endParameters;

    const auto& cov_1s_s = *(end_parameters_1s_s->covariance());

    // now it should be the same again, the PathScatter did the trick
    bool cov_s_approx = cov_1s_s.isApprox(cov_2s_s, 0.01);
    BOOST_TEST(cov_s_approx);
  }

}  // namespace Test
}  // namespace Acts
