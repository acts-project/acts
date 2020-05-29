// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <memory>
#include <random>

#include "Acts/EventData/FreeParameterSet.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

using namespace Acts::detail;

///
/// @brief Acts namespace
///
namespace Acts {
///
/// @brief Namespace for Acts unit tests
///
namespace Test {
/// @cond
namespace {
//~ // tolerance used for floating point comparison in this translation unit
//~ const double tol = 1e-6;

void check_residuals_for_bound_parameters() {
const double max = FreeParameterType<eFreeDir0>::max;
const double min = FreeParameterType<eFreeDir0>::min;
double tx_1 = 0.7;
double tx_2 = 0.4;
ActsVectorD<1> dTx;
dTx << (tx_1 - tx_2);

// both parameters inside bounds, difference is positive
FreeParameterSet<eFreeDir0> bound1(std::nullopt, tx_1);
FreeParameterSet<eFreeDir0> bound2(std::nullopt, tx_2);
CHECK_CLOSE_REL(bound1.residual(bound2), dTx, tol);

// both parameters inside bound, difference negative
dTx << (tx_2 - tx_1);
CHECK_CLOSE_REL(bound2.residual(bound1), dTheta, tol);

// one parameter above upper bound, difference positive
tx_1 = max + 1;
bound1.setParameter<eFreeDir0>(tx_1);
dTx << max - tx_2;
CHECK_CLOSE_REL(bound1.residual(bound2), dTx, tol);

// one parameter above upper bound, difference negative
dTx << tx_2 - max;
CHECK_CLOSE_REL(bound2.residual(bound1), dTx, tol);

// one parameter below lower bound, difference positive
tx_1 = min - 1;
bound1.setParameter<eFreeDir0>(tx_1);
dTx << tx_2 - min;
CHECK_CLOSE_REL(bound2.residual(bound1), dTx, tol);

// one parameter below lower bound, difference negative
dTx << min - tx_2;
CHECK_CLOSE_REL(bound1.residual(bound2), dTx, tol);

// both parameters outside bounds, both below
tx_1 = min - 1;
tx_2 = min - 2;
bound1.setParameter<eFreeDir0>(tx_1);
bound2.setParameter<eFreeDir0>(tx_2);
CHECK_SMALL(bound1.residual(bound2), tol);

// both parameters outside bounds, both above
tx_1 = max + 1;
tx_2 = max + 2;
bound1.setParameter<eFreeDir0>(tx_1);
bound2.setParameter<eFreeDir0>(tx_2);
CHECK_SMALL(bound1.residual(bound2), tol);

// both parameters outside bounds, one above, one below
tx_1 = max + 1;
tx_2 = min - 2;
bound1.setParameter<eFreeDir0>(tx_1);
bound2.setParameter<eFreeDir0>(tx_2);
dTx << max - min;
CHECK_CLOSE_REL(bound1.residual(bound2), dTx, tol);
dTx << min - max;
CHECK_CLOSE_REL(bound2.residual(bound1), dTx, tol);
}

void random_residual_tests() {
//~ // random number generators
//~ std::default_random_engine e;
//~ std::uniform_real_distribution<float> uniform_dist(-1000, 300);

//~ const double theta_max = BoundParameterType<eFreeDir2>::max;
//~ const double theta_min = BoundParameterType<eFreeDir2>::min;
//~ const double phi_max = BoundParameterType<eBoundPhi>::max;
//~ const double phi_min = BoundParameterType<eBoundPhi>::min;

//~ BoundVector parValues_1;
//~ BoundVector parValues_2;
//~ FullParameterSet parSet_1(std::nullopt, parValues_1);
//~ FullParameterSet parSet_2(std::nullopt, parValues_2);
//~ BoundVector residual;
//~ const unsigned int toys = 1000;
//~ for (unsigned int i = 0; i < toys; ++i) {
//~ const double loc0_1 = uniform_dist(e);
//~ const double loc1_1 = uniform_dist(e);
//~ const double phi_1 = uniform_dist(e);
//~ const double theta_1 = uniform_dist(e);
//~ const double qop_1 = uniform_dist(e);
//~ parValues_1 << loc0_1, loc1_1, phi_1, theta_1, qop_1, 0.;
//~ parSet_1.setParameters(parValues_1);

//~ const double loc0_2 = uniform_dist(e);
//~ const double loc1_2 = uniform_dist(e);
//~ const double phi_2 = uniform_dist(e);
//~ const double theta_2 = uniform_dist(e);
//~ const double qop_2 = uniform_dist(e);
//~ parValues_2 << loc0_2, loc1_2, phi_2, theta_2, qop_2, 0.;
//~ parSet_2.setParameters(parValues_2);

//~ const double delta_loc0 = loc0_1 - loc0_2;
//~ const double delta_loc1 = loc1_1 - loc1_2;
//~ // for theta make sure that the difference calculation considers the
//~ // restricted value range
//~ const double delta_theta =
//~ (theta_1 > theta_max ? theta_max
//~ : (theta_1 < theta_min ? theta_min : theta_1)) -
//~ (theta_2 > theta_max ? theta_max
//~ : (theta_2 < theta_min ? theta_min : theta_2));
//~ const double delta_qop = qop_1 - qop_2;
//~ residual = parSet_1.residual(parSet_2);

//~ // local parameters are unbound -> check for usual difference
//~ if (std::abs(residual(0) - delta_loc0) > tol) {
//~ BOOST_CHECK(false);
//~ break;
//~ }
//~ if (std::abs(residual(1) - delta_loc1) > tol) {
//~ BOOST_CHECK(false);
//~ break;
//~ }
//~ // phi is a cyclic parameter -> check that (unsigned) difference is not
//~ // larger than half period
//~ // check that corrected(corrected(phi_2) + residual) == corrected(phi_1)
//~ if (std::abs(get_cyclic_value(
//~ get_cyclic_value(phi_2, phi_min, phi_max) + residual(2),
//~ phi_min, phi_max) -
//~ get_cyclic_value(phi_1, phi_min, phi_max)) > tol or
//~ std::abs(residual(2)) > (phi_max - phi_min) / 2) {
//~ BOOST_CHECK(false);
//~ break;
//~ }
//~ // theta is bound -> check that (unsigned) difference is not larger then
//~ // allowed range, check corrected difference
//~ if (std::abs(residual(3) - delta_theta) > tol or
//~ std::abs(residual(3)) > (theta_max - theta_min)) {
//~ BOOST_CHECK(false);
//~ break;
//~ }
//~ // qop is unbound -> check usual difference
//~ if (std::abs(residual(4) - delta_qop) > tol) {
//~ BOOST_CHECK(false);
//~ break;
//~ }
//~ }
}
}  // namespace
/// @endcond

/**
 * @brief Unit test for checking consistency of FreeParameterSet class
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
BOOST_AUTO_TEST_CASE(free_parset_consistency_tests) {
  // check template parameter based information
  BOOST_CHECK((FreeParameterSet<eFreePos0, eFreePos1>::size() == 2));

  // covariance matrix
  ActsSymMatrixD<3> cov;
  cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;

  // parameter values
  double x = 0.5;
  double y = -0.2;
  double z = 0.3;
  ActsVectorD<3> parValues(x, y, z);

  // parameter set with covariance matrix
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> parSet_with_cov(cov, x, y,
                                                                    z);

  // check number and type of stored parameters
  BOOST_CHECK(parSet_with_cov.size() == 3);
  BOOST_CHECK(parSet_with_cov.contains<eFreePos0>());
  BOOST_CHECK(parSet_with_cov.contains<eFreePos1>());
  BOOST_CHECK(parSet_with_cov.contains<eFreePos2>());
  BOOST_CHECK(not parSet_with_cov.contains<eFreeTime>());
  BOOST_CHECK(not parSet_with_cov.contains<eFreeDir0>());
  BOOST_CHECK(not parSet_with_cov.contains<eFreeDir1>());
  BOOST_CHECK(not parSet_with_cov.contains<eFreeDir2>());
  BOOST_CHECK(not parSet_with_cov.contains<eFreeQOverP>());

  // check stored parameter values
  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos0>() == x);
  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos1>() == y);
  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos2>() == z);
  BOOST_CHECK(parSet_with_cov.getParameters() == parValues);

  // check stored covariance
  BOOST_CHECK(parSet_with_cov.getCovariance());
  BOOST_CHECK(*parSet_with_cov.getCovariance() == cov);
  BOOST_CHECK(parSet_with_cov.getUncertainty<eFreePos0>() == sqrt(cov(0, 0)));
  BOOST_CHECK(parSet_with_cov.getUncertainty<eFreePos1>() == sqrt(cov(1, 1)));
  BOOST_CHECK(parSet_with_cov.getUncertainty<eFreePos2>() == sqrt(cov(2, 2)));

  // same parameter set without covariance matrix
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> parSet_without_cov(
      std::nullopt, parValues);

  BOOST_CHECK(!parSet_without_cov.getCovariance());
  BOOST_CHECK(parSet_without_cov.getUncertainty<eFreePos0>() < 0);
  BOOST_CHECK(parSet_without_cov.getUncertainty<eFreePos1>() < 0);
  BOOST_CHECK(parSet_without_cov.getUncertainty<eFreePos2>() < 0);
  BOOST_CHECK(parSet_without_cov.getParameters() ==
              parSet_with_cov.getParameters());

  // set new covariance matrix
  parSet_without_cov.setCovariance(cov);

  BOOST_CHECK(parSet_without_cov.getCovariance());
  BOOST_CHECK(*parSet_without_cov.getCovariance() == cov);

  // set new parameter values
  double newX = 0.1;
  double newY = 0.6;
  double newZ = -0.15;
  parValues << newX, newY, newZ;
  parSet_with_cov.setParameter<eFreePos0>(newX);
  parSet_with_cov.setParameter<eFreePos1>(newY);
  parSet_with_cov.setParameter<eFreePos2>(newZ);

  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos0>() == newX);
  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos1>() == newY);
  BOOST_CHECK(parSet_with_cov.getParameter<eFreePos2>() == newZ);
  BOOST_CHECK(parSet_with_cov.getParameters() == parValues);
}

/**
 * @brief Unit test for copy/assignment/swap in FreeParameterSet class
 *
 * The behavior of the following functions is checked:
 * -# ParameterSet::ParameterSet
 * -# ParameterSet::operator=
 * -# ParameterSet::swap
 */
BOOST_AUTO_TEST_CASE(free_parset_copy_assignment_tests) {
  // covariance matrix
  ActsSymMatrixD<3> cov;
  cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;

  // parameter values
  double x = 0.5;
  double y = -0.2;
  double z = 0.3;
  ActsVectorD<3> first_parValues(x, y, z);

  // parameter set with covariance matrix
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> first(cov, x, y, z);

  // check copy constructor
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> copy(first);
  BOOST_CHECK(first == copy);

  // check move constructor
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> moved(std::move(copy));
  BOOST_CHECK(first == moved);

  // check assignment operator
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> assigned = moved;
  BOOST_CHECK(assigned == moved);

  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> other(std::nullopt, 0, 1.7,
                                                          -0.15);
  BOOST_CHECK(assigned != other);
  assigned = other;
  BOOST_CHECK(assigned == other);

  // check move assignment
  BOOST_CHECK(first != assigned);
  first = FreeParameterSet<eFreePos0, eFreePos1, eFreePos2>(assigned);
  BOOST_CHECK(first == assigned);

  // check swap method
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> lhs(cov, x, y, z);
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> rhs(std::nullopt, 2 * x,
                                                        2 * y, 2 * z);
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> lhs_copy = lhs;
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> rhs_copy = rhs;

  BOOST_CHECK(lhs != rhs && lhs == lhs_copy && rhs == rhs_copy);
  using std::swap;
  swap(lhs, rhs);
  BOOST_CHECK(lhs != rhs && rhs == lhs_copy && lhs == rhs_copy);
}

/**
 *  @brief Unit test for comparison operators in FreeParameterSet
 *
 *  @sa FreeParameterSet::operator==, FreeParameterSet::operator!=
 */
BOOST_AUTO_TEST_CASE(free_parset_comparison_tests) {
  // covariance matrix
  ActsSymMatrixD<3> cov;
  cov << 1, 0, 0, 0, 1.2, 0.2, 0, 0.2, 0.7;

  // parameter values
  double x = 0.5;
  double y = -0.2;
  double z = 0.3;

  // parameter set with covariance matrix
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> first(cov, x, y, z);
  FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> second(std::nullopt, 2 * x,
                                                           2 * y, 2 * z);

  // check self comparison
  BOOST_CHECK(first == first);
  BOOST_CHECK(not(first != first));

  // check mutual exclusivity
  BOOST_CHECK(first != second);
  BOOST_CHECK(not(first == second));
  first = second;
  BOOST_CHECK(first == second);

  // check that comparison fails for unequal parameter values
  second.setParameter<eFreePos0>(3 * x);
  BOOST_CHECK(first != second);
  first = second;
  BOOST_CHECK(first == second);

  second.setParameter<eFreePos1>(3 * y);
  BOOST_CHECK(first != second);
  first = second;
  BOOST_CHECK(first == second);

  second.setParameter<eFreePos2>(3 * z);
  BOOST_CHECK(first != second);
  first = second;
  BOOST_CHECK(first == second);

  // check that comparison fails for unequal covariance matrices
  second.setCovariance(cov);
  BOOST_CHECK(first != second);
  first = second;
  BOOST_CHECK(first == second);

  cov(0, 0) = 2 * cov(0, 0);
  second.setCovariance(cov);
  BOOST_CHECK(first != second);
  first = second;
  BOOST_CHECK(first == second);
}

/**
 * @brief Unit test for projection matrices in FreeParameterSet
 *
 * Checks the correctness of the projection matrices from the full parameter
 * space
 * onto different parameter sub-spaces
 *
 * @sa FreeParameterSet::projector
 */
BOOST_AUTO_TEST_CASE(free_parset_projection_tests) {
  ActsMatrixD<1, eFreeParametersSize> z_proj;
  z_proj << 0, 0, 1, 0, 0, 0, 0, 0;

  ActsMatrixD<2, eFreeParametersSize> x_qop_proj;
  x_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  ActsMatrixD<2, eFreeParametersSize> y_tz_proj;
  y_tz_proj << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;

  ActsMatrixD<3, eFreeParametersSize> x_y_z_proj;
  x_y_z_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0;

  ActsMatrixD<4, eFreeParametersSize> x_z_tz_qop_proj;
  x_z_tz_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  ActsMatrixD<5, eFreeParametersSize> x_y_z_tz_qop_proj;
  x_y_z_tz_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  ActsMatrixD<6, eFreeParametersSize> x_y_z_t_tz_qop_proj;
  x_y_z_t_tz_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1;

  ActsMatrixD<7, eFreeParametersSize> x_y_z_t_ty_tz_qop_proj;
  x_y_z_t_ty_tz_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  ActsMatrixD<eFreeParametersSize, eFreeParametersSize>
      x_y_z_t_tx_ty_tz_qop_proj;
  x_y_z_t_tx_ty_tz_qop_proj << 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1;

  BOOST_CHECK((FreeParameterSet<eFreePos2>::projector() == z_proj));
  BOOST_CHECK(
      (FreeParameterSet<eFreePos0, eFreeQOverP>::projector() == x_qop_proj));
  BOOST_CHECK(
      (FreeParameterSet<eFreePos1, eFreeDir2>::projector() == y_tz_proj));
  BOOST_CHECK((FreeParameterSet<eFreePos0, eFreePos1, eFreePos2>::projector() ==
               x_y_z_proj));
  BOOST_CHECK((FreeParameterSet<eFreePos0, eFreePos2, eFreeDir2,
                                eFreeQOverP>::projector() == x_z_tz_qop_proj));
  BOOST_CHECK(
      (FreeParameterSet<eFreePos0, eFreePos1, eFreePos2, eFreeDir2,
                        eFreeQOverP>::projector() == x_y_z_tz_qop_proj));
  BOOST_CHECK(
      (FreeParameterSet<eFreePos0, eFreePos1, eFreePos2, eFreeTime, eFreeDir2, eFreeQOverP
                        >::projector() == x_y_z_t_tz_qop_proj));
  BOOST_CHECK((FreeParameterSet<eFreePos0, eFreePos1, eFreePos2, eFreeTime, eFreeDir1,
                                eFreeDir2, eFreeQOverP>::projector() ==
               x_y_z_t_ty_tz_qop_proj));
  BOOST_CHECK(
      (FreeParameterSet<eFreePos0, eFreePos1, eFreePos2, eFreeTime, eFreeDir0, eFreeDir1,
                        eFreeDir2, eFreeQOverP>::projector() ==
       x_y_z_t_tx_ty_tz_qop_proj));
}

/**
 * @brief Unit test for residuals between different FreeParameterSet objects
 *
 * The result of the residual calculation between two FreeParameterSet objects is
 * checked.
 * A test of the automatic correction of stored parameter values for
 * out-of-bounds
 * corrections is also implemented.
 *
* @sa FreeParameterSet::residual, FreeParameterSet::getParameter
*/
BOOST_AUTO_TEST_CASE(free_parset_residual_tests) {
// check unbound parameter type
const double large_number = 12443534120;
const double small_number = -924342675;
const double normal_number = 0.1234;
FreeParameterSet<eFreePos0, eFreePos1, eFreeQOverP> unbound(
std::nullopt, small_number, large_number, normal_number);
BOOST_CHECK(unbound.getParameter<eFreePos0>() == small_number);
BOOST_CHECK(unbound.getParameter<eFreePos1>() == large_number);
BOOST_CHECK(unbound.getParameter<eFreeQOverP>() == normal_number);

// check bound parameter type
FreeParameterSet<eFreeDir2> bound(std::nullopt, small_number);
BOOST_CHECK((bound.getParameter<eFreeDir2>() ==
FreeParameterType<eFreeDir2>::min));
bound.setParameter<eFreeDir2>(large_number);
BOOST_CHECK((bound.getParameter<eFreeDir2>() ==
FreeParameterType<eFreeDir2>::max));
bound.setParameter<eFreeDir2>(normal_number);
BOOST_CHECK((bound.getParameter<eFreeDir2>() == normal_number));

// check residual calculation

// input numbers
const double first_x = 0.3;
const double first_y = 0.9;
const double first_z = 0.7;

const double second_x = 2.7;
const double second_y = -0.9;
const double second_z = 0.35;

// expected results for residual second wrt first
const double delta_x = second_x - first_x;
const double delta_y = second_y - first_y;
const double delta_z = second_z - first_z;
ActsVectorD<3> residuals(delta_x, delta_y, delta_z);

FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> first(
std::nullopt, first_x, first_y, first_z);
FreeParameterSet<eFreePos0, eFreePos1, eFreePos2> second(
std::nullopt, second_x, second_y, second_z);
CHECK_CLOSE_REL(residuals, second.residual(first), 1e-6);

// some more checks for bound variables
check_residuals_for_bound_parameters();

// // inspecific residual tests with random numbers
// random_residual_tests();
}

//~ template <ParID_t... params>
//~ using ParSet = FreeParameterSet<params...>;

//~ /**
//~ * @brief Unit test for index-/type-based access of coordinates
//~ *
//~ * @sa FreeParameterSet::getIndex
//~ * @sa FreeParameterSet::getParID
//~ */
//~ BOOST_AUTO_TEST_CASE(parset_parID_mapping) {
//~ // check logic for type-based access
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getIndex<eFreePos0>() == 0));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getIndex<eFreePos1>() == 1));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getIndex<eBoundPhi>() == 2));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getIndex<eBoundQOverP>() == 3));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getIndex<eT>() == 4));

//~ // check logic for index-based access
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getParID<0>() == eFreePos0));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getParID<1>() == eFreePos1));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getParID<2>() == eBoundPhi));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getParID<3>() == eBoundQOverP));
//~ BOOST_CHECK((ParSet<eFreePos0, eFreePos1, eBoundPhi, eBoundQOverP,
//~ eT>::getParID<4>() == eT));

//~ // check consistency
//~ using FullSet = FullFreeParameterSet;
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<0>()>() == 0));
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<1>()>() == 1));
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<2>()>() == 2));
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<3>()>() == 3));
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<4>()>() == 4));
//~ BOOST_CHECK((FullSet::getIndex<FullSet::getParID<5>()>() == 5));

//~ BOOST_CHECK(
//~ (FullSet::getParID<FullSet::getIndex<eFreePos0>()>() == eFreePos0));
//~ BOOST_CHECK(
//~ (FullSet::getParID<FullSet::getIndex<eFreePos1>()>() == eFreePos1));
//~ BOOST_CHECK(
//~ (FullSet::getParID<FullSet::getIndex<eBoundPhi>()>() == eBoundPhi));
//~ BOOST_CHECK(
//~ (FullSet::getParID<FullSet::getIndex<eFreeDir2>()>() == eFreeDir2));
//~ BOOST_CHECK(
//~ (FullSet::getParID<FullSet::getIndex<eBoundQOverP>()>() == eBoundQOverP));
//~ BOOST_CHECK((FullSet::getParID<FullSet::getIndex<eT>()>() == eT));

//~ // consistency of types
//~ BOOST_CHECK((std::is_same<std::remove_cv<decltype(
//~ at_index<ParID_t, 0, eFreePos0>::value)>::type,
//~ decltype(eFreePos0)>::value));
//~ BOOST_CHECK((std::is_same<decltype(FullSet::getParID<0>()),
//~ decltype(eFreePos0)>::value));
//~ }
}  // namespace Test
}  // namespace Acts
