// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE ParameterSet Tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include <cmath>
#include <memory>
#include <random>

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

using namespace Acts::detail;

///
/// @brief Ats namespace
///
namespace Acts {
///
/// @brief Namespace for ATS unit tests
///
namespace Test {
  /// @cond
  namespace {
    // tolerance used for floating point comparison in this translation unit
    const double tol = 1e-6;

    double
    get_cyclic_value(double value, double min, double max)
    {
      return value - (max - min) * std::floor((value - min) / (max - min));
    }

    double
    get_cyclic_difference(double a, double b, double min, double max)
    {
      const double period      = max - min;
      const double half_period = period / 2;
      a                        = get_cyclic_value(a, min, max);
      b                        = get_cyclic_value(b, min, max);
      double raw_diff          = a - b;
      double diff              = (raw_diff > half_period)
          ? raw_diff - period
          : ((raw_diff < -half_period) ? period + raw_diff : raw_diff);
      return diff;
    }

    void
    check_residuals_for_bound_parameters()
    {
      const double   max     = par_type_t<ParID_t::eTHETA>::max;
      const double   min     = par_type_t<ParID_t::eTHETA>::min;
      double         theta_1 = 0.7 * M_PI;
      double         theta_2 = 0.4 * M_PI;
      ActsVectorD<1> dTheta;
      dTheta << (theta_1 - theta_2);

      // both parameters inside bounds, difference is positive
      ParameterSet<ParID_t::eTHETA> bound1(nullptr, theta_1);
      ParameterSet<ParID_t::eTHETA> bound2(nullptr, theta_2);
      CHECK_CLOSE_REL(bound1.residual(bound2), dTheta, tol);

      // both parameters inside bound, difference negative
      dTheta << (theta_2 - theta_1);
      CHECK_CLOSE_REL(bound2.residual(bound1), dTheta, tol);

      // one parameter above upper bound, difference positive
      theta_1 = max + 1;
      bound1.setParameter<ParID_t::eTHETA>(theta_1);
      dTheta << max - theta_2;
      CHECK_CLOSE_REL(bound1.residual(bound2), dTheta, tol);

      // one parameter above upper bound, difference negative
      dTheta << theta_2 - max;
      CHECK_CLOSE_REL(bound2.residual(bound1), dTheta, tol);

      // one parameter below lower bound, difference positive
      theta_1 = min - 1;
      bound1.setParameter<ParID_t::eTHETA>(theta_1);
      dTheta << theta_2 - min;
      CHECK_CLOSE_REL(bound2.residual(bound1), dTheta, tol);

      // one parameter below lower bound, difference negative
      dTheta << min - theta_2;
      CHECK_CLOSE_REL(bound1.residual(bound2), dTheta, tol);

      // both parameters outside bounds, both below
      theta_1 = min - 1;
      theta_2 = min - 2;
      bound1.setParameter<ParID_t::eTHETA>(theta_1);
      bound2.setParameter<ParID_t::eTHETA>(theta_2);
      CHECK_SMALL(bound1.residual(bound2), tol);

      // both parameters outside bounds, both above
      theta_1 = max + 1;
      theta_2 = max + 2;
      bound1.setParameter<ParID_t::eTHETA>(theta_1);
      bound2.setParameter<ParID_t::eTHETA>(theta_2);
      CHECK_SMALL(bound1.residual(bound2), tol);

      // both parameters outside bounds, one above, one below
      theta_1 = max + 1;
      theta_2 = min - 2;
      bound1.setParameter<ParID_t::eTHETA>(theta_1);
      bound2.setParameter<ParID_t::eTHETA>(theta_2);
      dTheta << max - min;
      CHECK_CLOSE_REL(bound1.residual(bound2), dTheta, tol);
      dTheta << min - max;
      CHECK_CLOSE_REL(bound2.residual(bound1), dTheta, tol);
    }

    void
    check_residuals_for_cyclic_parameters()
    {
      const double max = par_type_t<ParID_t::ePHI>::max;
      const double min = par_type_t<ParID_t::ePHI>::min;

      double         phi_1 = 0.7 * M_PI;
      double         phi_2 = 0.4 * M_PI;
      ActsVectorD<1> dPhi;
      dPhi << (phi_1 - phi_2);

      ParameterSet<ParID_t::ePHI> cyclic1(nullptr, phi_1);
      ParameterSet<ParID_t::ePHI> cyclic2(nullptr, phi_2);

      // no boundary crossing, difference is positive
      CHECK_CLOSE_REL(cyclic1.residual(cyclic2), dPhi, tol);

      // no boundary crossing, difference is negative
      CHECK_CLOSE_REL(cyclic2.residual(cyclic1), -dPhi, tol);

      // forward boundary crossing
      phi_1 = -0.9 * M_PI;
      cyclic1.setParameter<ParID_t::ePHI>(phi_1);
      dPhi << get_cyclic_difference(phi_1, phi_2, min, max);
      CHECK_CLOSE_REL(cyclic1.residual(cyclic2), dPhi, tol);
      CHECK_CLOSE_REL(cyclic2.residual(cyclic1), -dPhi, tol);

      // backward boundary crossing
      phi_1 = 0.7 * M_PI;
      phi_2 = -0.9 * M_PI;
      cyclic1.setParameter<ParID_t::ePHI>(phi_1);
      cyclic2.setParameter<ParID_t::ePHI>(phi_2);
      dPhi << get_cyclic_difference(phi_1, phi_2, min, max);
      CHECK_CLOSE_REL(cyclic1.residual(cyclic2), dPhi, tol);
      CHECK_CLOSE_REL(cyclic2.residual(cyclic1), -dPhi, tol);
    }

    void
    random_residual_tests()
    {
      // random number generators
      std::default_random_engine            e;
      std::uniform_real_distribution<float> uniform_dist(-1000, 300);

      const double theta_max = par_type_t<ParID_t::eTHETA>::max;
      const double theta_min = par_type_t<ParID_t::eTHETA>::min;
      const double phi_max   = par_type_t<ParID_t::ePHI>::max;
      const double phi_min   = par_type_t<ParID_t::ePHI>::min;

      ActsVectorD<5>     parValues_1;
      ActsVectorD<5>     parValues_2;
      FullParameterSet   parSet_1(nullptr, parValues_1);
      FullParameterSet   parSet_2(nullptr, parValues_2);
      ActsVectorD<5>     residual;
      const unsigned int toys = 1000;
      for (unsigned int i = 0; i < toys; ++i) {
        const double loc0_1  = uniform_dist(e);
        const double loc1_1  = uniform_dist(e);
        const double phi_1   = uniform_dist(e);
        const double theta_1 = uniform_dist(e);
        const double qop_1   = uniform_dist(e);
        parValues_1 << loc0_1, loc1_1, phi_1, theta_1, qop_1;
        parSet_1.setParameters(parValues_1);

        const double loc0_2  = uniform_dist(e);
        const double loc1_2  = uniform_dist(e);
        const double phi_2   = uniform_dist(e);
        const double theta_2 = uniform_dist(e);
        const double qop_2   = uniform_dist(e);
        parValues_2 << loc0_2, loc1_2, phi_2, theta_2, qop_2;
        parSet_2.setParameters(parValues_2);

        const double delta_loc0 = loc0_1 - loc0_2;
        const double delta_loc1 = loc1_1 - loc1_2;
        // for theta make sure that the difference calculation considers the
        // restricted value range
        const double delta_theta
            = (theta_1 > theta_max
                   ? theta_max
                   : (theta_1 < theta_min ? theta_min : theta_1))
            - (theta_2 > theta_max
                   ? theta_max
                   : (theta_2 < theta_min ? theta_min : theta_2));
        const double delta_qop = qop_1 - qop_2;
        residual               = parSet_1.residual(parSet_2);

        // local parameters are unbound -> check for usual difference
        if (std::abs(residual(0) - delta_loc0) > tol) {
          BOOST_CHECK(false);
          break;
        }
        if (std::abs(residual(1) - delta_loc1) > tol) {
          BOOST_CHECK(false);
          break;
        }
        // phi is a cyclic parameter -> check that (unsigned) difference is not
        // larger than half period
        // check that corrected(corrected(phi_2) + residual) == corrected(phi_1)
        if (std::abs(get_cyclic_value(get_cyclic_value(phi_2, phi_min, phi_max)
                                          + residual(2),
                                      phi_min,
                                      phi_max)
                     - get_cyclic_value(phi_1, phi_min, phi_max))
                > tol
            or std::abs(residual(2)) > (phi_max - phi_min) / 2) {
          BOOST_CHECK(false);
          break;
        }
        // theta is bound -> check that (unsigned) difference is not larger then
        // allowed range, check corrected difference
        if (std::abs(residual(3) - delta_theta) > tol
            or std::abs(residual(3)) > (theta_max - theta_min)) {
          BOOST_CHECK(false);
          break;
        }
        // qop is unbound -> check usual difference
        if (std::abs(residual(4) - delta_qop) > tol) {
          BOOST_CHECK(false);
          break;
        }
      }
    }
  }
  /// @endcond

  /**
   * @brief Unit test for checking consistency of ParameterSet class
   *
   * The following functions are tested to yield the expected result/behavior:
   * -# ParameterSet::size
   * -# ParameterSet::contains
   * -# ParameterSet::getParameter
   * -# ParameterSet::getParameters
   * -# ParameterSet::getCovariance
   * -# ParameterSet::setCovariance
   * -# ParameterSet::setParameter
   * -# ParameterSet::getUncertainty
   */
  BOOST_AUTO_TEST_CASE(parset_consistency_tests)
  {
    // check template parameter based information
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1>::size() == 2));

    // covariance matrix
    ActsSymMatrixD<3> cov;
    cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;
    auto pCovMatrix = std::make_unique<const ActsSymMatrixD<3>>(cov);

    // parameter values
    double loc0 = 0.5;
    double loc1 = -0.2;
    double phi  = 0.3 * M_PI;  // this should be within [-M_PI,M_PI) to avoid
                               // failed tests due to angle range corrections
    ActsVectorD<3> parValues(loc0, loc1, phi);

    // parameter set with covariance matrix
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI>
        parSet_with_cov(std::move(pCovMatrix), loc0, loc1, phi);

    // check number and type of stored parameters
    BOOST_CHECK(parSet_with_cov.size() == 3);
    BOOST_CHECK(parSet_with_cov.contains<ParID_t::eLOC_0>());
    BOOST_CHECK(parSet_with_cov.contains<ParID_t::eLOC_1>());
    BOOST_CHECK(parSet_with_cov.contains<ParID_t::ePHI>());
    BOOST_CHECK(not parSet_with_cov.contains<ParID_t::eTHETA>());
    BOOST_CHECK(not parSet_with_cov.contains<ParID_t::eQOP>());

    // check stored parameter values
    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::eLOC_0>() == loc0);
    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::eLOC_1>() == loc1);
    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::ePHI>() == phi);
    BOOST_CHECK(parSet_with_cov.getParameters() == parValues);

    // check stored covariance
    BOOST_CHECK(parSet_with_cov.getCovariance() != 0);
    BOOST_CHECK(*parSet_with_cov.getCovariance() == cov);
    BOOST_CHECK(parSet_with_cov.getUncertainty<ParID_t::eLOC_0>()
                == sqrt(cov(0, 0)));
    BOOST_CHECK(parSet_with_cov.getUncertainty<ParID_t::eLOC_1>()
                == sqrt(cov(1, 1)));
    BOOST_CHECK(parSet_with_cov.getUncertainty<ParID_t::ePHI>()
                == sqrt(cov(2, 2)));

    // same parameter set without covariance matrix
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI>
        parSet_without_cov(nullptr, parValues);

    BOOST_CHECK(parSet_without_cov.getCovariance() == 0);
    BOOST_CHECK(parSet_without_cov.getUncertainty<ParID_t::eLOC_0>() < 0);
    BOOST_CHECK(parSet_without_cov.getUncertainty<ParID_t::eLOC_1>() < 0);
    BOOST_CHECK(parSet_without_cov.getUncertainty<ParID_t::ePHI>() < 0);
    BOOST_CHECK(parSet_without_cov.getParameters()
                == parSet_with_cov.getParameters());

    // set new covariance matrix
    parSet_without_cov.setCovariance(
        std::make_unique<const ActsSymMatrixD<3>>(cov));

    BOOST_CHECK(parSet_without_cov.getCovariance() != 0);
    BOOST_CHECK(*parSet_without_cov.getCovariance() == cov);

    // set new parameter values
    double newLoc0 = 0.1;
    double newLoc1 = 0.6;
    double newPhi  = -0.15 * M_PI;
    parValues << newLoc0, newLoc1, newPhi;
    parSet_with_cov.setParameter<ParID_t::eLOC_0>(newLoc0);
    parSet_with_cov.setParameter<ParID_t::eLOC_1>(newLoc1);
    parSet_with_cov.setParameter<ParID_t::ePHI>(newPhi);

    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::eLOC_0>() == newLoc0);
    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::eLOC_1>() == newLoc1);
    BOOST_CHECK(parSet_with_cov.getParameter<ParID_t::ePHI>() == newPhi);
    BOOST_CHECK(parSet_with_cov.getParameters() == parValues);
  }

  /**
   * @brief Unit test for copy/assignment/swap in ParameterSet class
   *
   * The behavior of the following functions is checked:
   * -# ParameterSet::ParameterSet
   * -# ParameterSet::operator=
   * -# ParameterSet::swap
   */
  BOOST_AUTO_TEST_CASE(parset_copy_assignment_tests)
  {
    // covariance matrix
    ActsSymMatrixD<3> cov;
    cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;
    auto pCovMatrix = std::make_unique<const ActsSymMatrixD<3>>(cov);

    // parameter values
    double loc0 = 0.5;
    double loc1 = -0.2;
    double phi  = 0.3 * M_PI;  // this should be within [-M_PI,M_PI) to avoid
                               // failed tests due to angle range corrections
    ActsVectorD<3> first_parValues(loc0, loc1, phi);

    // parameter set with covariance matrix
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> first(
        std::move(pCovMatrix), loc0, loc1, phi);

    // check copy constructor
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> copy(first);
    BOOST_CHECK(first == copy);

    // check move constructor
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> moved(
        std::move(copy));
    BOOST_CHECK(first == moved);

    // check assignment operator
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> assigned
        = moved;
    BOOST_CHECK(assigned == moved);

    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> other(
        nullptr, 0, 1.7, -0.15);
    BOOST_CHECK(assigned != other);
    assigned = other;
    BOOST_CHECK(assigned == other);

    // check move assignment
    BOOST_CHECK(first != assigned);
    first = ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI>(
        assigned);
    BOOST_CHECK(first == assigned);

    // check swap method
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> lhs(
        std::move(pCovMatrix), loc0, loc1, phi);
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> rhs(
        nullptr, 2 * loc0, 2 * loc1, 2 * phi);
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> lhs_copy
        = lhs;
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> rhs_copy
        = rhs;

    BOOST_CHECK(lhs != rhs && lhs == lhs_copy && rhs == rhs_copy);
    using std::swap;
    swap(lhs, rhs);
    BOOST_CHECK(lhs != rhs && rhs == lhs_copy && lhs == rhs_copy);
  }

  /**
   *  @brief Unit test for comparison operators in ParameterSet
   *
   *  @sa ParameterSet::operator==, ParameterSet::operator!=
   */
  BOOST_AUTO_TEST_CASE(parset_comparison_tests)
  {
    // covariance matrix
    ActsSymMatrixD<3> cov;
    cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;
    auto pCovMatrix = std::make_unique<const ActsSymMatrixD<3>>(cov);

    // parameter values
    double loc0 = 0.5;
    double loc1 = -0.2;
    double phi  = 0.3 * M_PI;  // this should be within [-M_PI,M_PI) to avoid
                               // failed tests due to angle range corrections

    // parameter set with covariance matrix
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> first(
        std::move(pCovMatrix), loc0, loc1, phi);
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI> second(
        nullptr, 2 * loc0, 2 * loc1, 2 * phi);

    // check self comparison
    BOOST_CHECK(first == first);
    BOOST_CHECK(not(first != first));

    // check mutual exclusivity
    BOOST_CHECK(first != second);
    BOOST_CHECK(not(first == second));
    first = second;
    BOOST_CHECK(first == second);

    // check that comparison fails for unequal parameter values
    second.setParameter<ParID_t::eLOC_0>(3 * loc0);
    BOOST_CHECK(first != second);
    first = second;
    BOOST_CHECK(first == second);

    second.setParameter<ParID_t::eLOC_1>(3 * loc1);
    BOOST_CHECK(first != second);
    first = second;
    BOOST_CHECK(first == second);

    second.setParameter<ParID_t::ePHI>(3 * phi);
    BOOST_CHECK(first != second);
    first = second;
    BOOST_CHECK(first == second);

    // check that comparison fails for unequal covariance matrices
    second.setCovariance(std::make_unique<const ActsSymMatrixD<3>>(cov));
    BOOST_CHECK(first != second);
    first = second;
    BOOST_CHECK(first == second);

    cov(0, 0) = 2 * cov(0, 0);
    second.setCovariance(std::make_unique<const ActsSymMatrixD<3>>(cov));
    BOOST_CHECK(first != second);
    first = second;
    BOOST_CHECK(first == second);
  }

  /**
   * @brief Unit test for projection matrices in ParameterSet
   *
   * Checks the correctness of the projection matrices from the full parameter
   * space
   * onto different parameter sub-spaces
   *
   * @sa ParameterSet::projector
   */
  BOOST_AUTO_TEST_CASE(parset_projection_tests)
  {
    ActsMatrixD<1, 5> phi_proj;
    phi_proj << 0, 0, 1, 0, 0;

    ActsMatrixD<2, 5> loc0_qop_proj;
    loc0_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

    ActsMatrixD<2, 5> loc1_theta_proj;
    loc1_theta_proj << 0, 1, 0, 0, 0, 0, 0, 0, 1, 0;

    ActsMatrixD<3, 5> loc0_loc1_phi_proj;
    loc0_loc1_phi_proj << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0;

    ActsMatrixD<4, 5> loc0_phi_theta_qop_proj;
    loc0_phi_theta_qop_proj << 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 1;

    ActsMatrixD<5, 5> loc0_loc1_phi_theta_qop_proj;
    loc0_loc1_phi_theta_qop_proj << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 1;

    BOOST_CHECK((ParameterSet<ParID_t::ePHI>::projector() == phi_proj));
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_0, ParID_t::eQOP>::projector()
                 == loc0_qop_proj));
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_1, ParID_t::eTHETA>::projector()
                 == loc1_theta_proj));
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::ePHI>::
                     projector()
                 == loc0_loc1_phi_proj));
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_0,
                              ParID_t::ePHI,
                              ParID_t::eTHETA,
                              ParID_t::eQOP>::projector()
                 == loc0_phi_theta_qop_proj));
    BOOST_CHECK((ParameterSet<ParID_t::eLOC_0,
                              ParID_t::eLOC_1,
                              ParID_t::ePHI,
                              ParID_t::eTHETA,
                              ParID_t::eQOP>::projector()
                 == loc0_loc1_phi_theta_qop_proj));
  }

  /**
   * @brief Unit test for residuals between different ParameterSet objects
   *
   * The result of the residual calculation between two ParameterSet objects is
   * checked.
   * A test of the automatic correction of stored parameter values for
   * out-of-bounds/cyclic
   * corrections is also implemented.
   *
   * @sa ParameterSet::residual, ParameterSet::getParameter
   */
  BOOST_AUTO_TEST_CASE(parset_residual_tests)
  {
    // check unbound parameter type
    const double large_number  = 12443534120;
    const double small_number  = -924342675;
    const double normal_number = 1.234;
    ParameterSet<ParID_t::eLOC_0, ParID_t::eLOC_1, ParID_t::eQOP> unbound(
        nullptr, small_number, large_number, normal_number);
    BOOST_CHECK(unbound.getParameter<ParID_t::eLOC_0>() == small_number);
    BOOST_CHECK(unbound.getParameter<ParID_t::eLOC_1>() == large_number);
    BOOST_CHECK(unbound.getParameter<ParID_t::eQOP>() == normal_number);

    // check bound parameter type
    ParameterSet<ParID_t::eTHETA> bound(nullptr, small_number);
    BOOST_CHECK((bound.getParameter<ParID_t::eTHETA>()
                 == par_type_t<ParID_t::eTHETA>::min));
    bound.setParameter<ParID_t::eTHETA>(large_number);
    BOOST_CHECK((bound.getParameter<ParID_t::eTHETA>()
                 == par_type_t<ParID_t::eTHETA>::max));
    bound.setParameter<ParID_t::eTHETA>(normal_number);
    BOOST_CHECK((bound.getParameter<ParID_t::eTHETA>() == normal_number));

    // check cyclic parameter type
    ParameterSet<ParID_t::ePHI> cyclic(nullptr, small_number);
    // calculate expected results
    const double min = par_type_t<ParID_t::ePHI>::min;
    const double max = par_type_t<ParID_t::ePHI>::max;
    // check that difference between original phi and stored phi is a multiple
    // of the cyclic period
    double multiple
        = (cyclic.getParameter<ParID_t::ePHI>() - small_number) / (max - min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() >= min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() < max);
    BOOST_CHECK(std::abs(multiple - std::floor(multiple + 0.5)) < tol);

    cyclic.setParameter<ParID_t::ePHI>(large_number);
    multiple
        = (cyclic.getParameter<ParID_t::ePHI>() - large_number) / (max - min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() >= min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() < max);
    BOOST_CHECK(std::abs(multiple - std::floor(multiple + 0.5)) < tol);

    cyclic.setParameter<ParID_t::ePHI>(normal_number);
    multiple
        = (cyclic.getParameter<ParID_t::ePHI>() - normal_number) / (max - min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() >= min);
    BOOST_CHECK(cyclic.getParameter<ParID_t::ePHI>() < max);
    BOOST_CHECK(std::abs(multiple - std::floor(multiple + 0.5)) < tol);

    // check residual calculation

    // input numbers
    const double first_loc0  = 0.3;
    const double first_phi   = 0.9 * M_PI;
    const double first_theta = 0.7 * M_PI;

    const double second_loc0  = 2.7;
    const double second_phi   = -0.9 * M_PI;
    const double second_theta = 0.35 * M_PI;

    // expected results for residual second wrt first
    const double delta_loc0 = second_loc0 - first_loc0;
    const double delta_phi
        = get_cyclic_difference(second_phi, first_phi, min, max);
    const double   delta_theta = second_theta - first_theta;
    ActsVectorD<3> residuals(delta_loc0, delta_phi, delta_theta);

    ParameterSet<ParID_t::eLOC_0, ParID_t::ePHI, ParID_t::eTHETA> first(
        nullptr, first_loc0, first_phi, first_theta);
    ParameterSet<ParID_t::eLOC_0, ParID_t::ePHI, ParID_t::eTHETA> second(
        nullptr, second_loc0, second_phi, second_theta);
    CHECK_CLOSE_REL(residuals, second.residual(first), 1e-6);

    // some more checks for bound variables
    check_residuals_for_bound_parameters();

    // some more checks for cyclic variables
    check_residuals_for_cyclic_parameters();

    // inspecific residual tests with random numbers
    random_residual_tests();
  }

  template <ParID_t... params>
  using ParSet = ParameterSet<params...>;

  /**
   * @brief Unit test for index-/type-based access of coordinates
   *
   * @sa ParameterSet::getIndex
   * @sa ParameterSet::getParID
   */
  BOOST_AUTO_TEST_CASE(parset_parID_mapping)
  {
    // check logic for type-based access
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getIndex<eLOC_0>() == 0));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getIndex<eLOC_1>() == 1));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getIndex<ePHI>() == 2));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getIndex<eQOP>() == 3));

    // check logic for index-based access
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getParID<0>() == eLOC_0));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getParID<1>() == eLOC_1));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getParID<2>() == ePHI));
    BOOST_CHECK((ParSet<eLOC_0, eLOC_1, ePHI, eQOP>::getParID<3>() == eQOP));

    // check consistency
    using FullSet = FullParameterSet;
    BOOST_CHECK((FullSet::getIndex<FullSet::getParID<0>()>() == 0));
    BOOST_CHECK((FullSet::getIndex<FullSet::getParID<1>()>() == 1));
    BOOST_CHECK((FullSet::getIndex<FullSet::getParID<2>()>() == 2));
    BOOST_CHECK((FullSet::getIndex<FullSet::getParID<3>()>() == 3));
    BOOST_CHECK((FullSet::getIndex<FullSet::getParID<4>()>() == 4));

    BOOST_CHECK((FullSet::getParID<FullSet::getIndex<eLOC_0>()>() == eLOC_0));
    BOOST_CHECK((FullSet::getParID<FullSet::getIndex<eLOC_1>()>() == eLOC_1));
    BOOST_CHECK((FullSet::getParID<FullSet::getIndex<ePHI>()>() == ePHI));
    BOOST_CHECK((FullSet::getParID<FullSet::getIndex<eTHETA>()>() == eTHETA));
    BOOST_CHECK((FullSet::getParID<FullSet::getIndex<eQOP>()>() == eQOP));

    // consistency of types
    BOOST_CHECK((std::is_same<std::remove_cv<decltype(
                                  at_index<ParID_t, 0, eLOC_0>::value)>::type,
                              decltype(eLOC_0)>::value));
    BOOST_CHECK((std::is_same<decltype(FullSet::getParID<0>()),
                              decltype(eLOC_0)>::value));
  }
}  // namespace Test
}  // namespace Acts
