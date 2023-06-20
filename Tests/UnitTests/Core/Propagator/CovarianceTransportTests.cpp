// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <memory>
#include <optional>
#include <tuple>
#include <variant>

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_CASE(covariance_transport_invalid) {
  CovarianceCache covCache;
  BOOST_CHECK(covCache.applyTransport == false);
}

BOOST_AUTO_TEST_CASE(covariance_transport_bound_start) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  Vector3 position{1., 2., 3.};
  ActsScalar time = 4.;
  ActsScalar qop = 0.125;
  Vector3 direction = Vector3(5., 6., 7.).normalized();
  auto planeSurface = Surface::makeShared<PlaneSurface>(position, direction);
  auto otherSurface = Surface::makeShared<PlaneSurface>(
      position, Vector3(6., 7., 8.).normalized());

  // Free & bound parameters
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;

  BoundVector boundParameters;
  boundParameters << 0., 0., VectorHelpers::phi(direction),
      VectorHelpers::theta(direction), qop, time;

  BoundSymMatrix boundCovariance = 8. * BoundSymMatrix::Identity();

  CovarianceCache covCache(tgContext, *planeSurface, position, boundParameters,
                           boundCovariance);
  // Test that the transport should be applied now
  BOOST_CHECK(covCache.applyTransport == true);
  // Test that the set covariance is 5x5 and what it has been set
  BOOST_CHECK_EQUAL(std::get<BoundSymMatrix>(covCache.covariance),
                    boundCovariance);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(covCache.covariance),
                    std::bad_variant_access);
  // Test that the bound to free jacobian has a value
  BOOST_CHECK(covCache.boundToFreeJacobian.has_value());
  BOOST_CHECK(not covCache.boundToFreeJacobian.value().isApprox(
      BoundToFreeMatrix::Zero()));
  // Test that the free transport jacobian, derivative is properly set up
  BOOST_CHECK_EQUAL(covCache.freeTransportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(covCache.freeToPathDerivatives, FreeVector::Zero());
  // Check that the surface is set
  BOOST_CHECK_EQUAL(covCache.atSurface, planeSurface.get());
  BOOST_CHECK_EQUAL(covCache.atPosition, position);

  // Transport bound to the same bound surface
  auto covJacAtBound = transportCovarianceToBound(tgContext, *planeSurface,
                                                  freeParameters, covCache);
  BOOST_CHECK(
      boundCovariance.isApprox(std::get<BoundSymMatrix>(covCache.covariance)));
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(covCache.covariance),
                    std::bad_variant_access);

  auto variantCovariance = std::get<0>(covJacAtBound);
  auto variantJacobian = std::get<1>(covJacAtBound);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to curvilinear on the same place
  auto covJacAtCurv = transportCovarianceToCurvilinear(direction, covCache);
  BOOST_CHECK(
      boundCovariance.isApprox(std::get<BoundSymMatrix>(covCache.covariance)));
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(covCache.covariance),
                    std::bad_variant_access);

  variantCovariance = std::get<0>(covJacAtCurv);
  variantJacobian = std::get<1>(covJacAtCurv);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to Free on the same place
  auto covJacAtFree = transportCovarianceToFree(covCache);
  variantCovariance = std::get<0>(covJacAtFree);
  variantJacobian = std::get<1>(covJacAtFree);
  BOOST_CHECK_THROW(std::get<BoundSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to another bound surface
  covJacAtBound = transportCovarianceToBound(tgContext, *otherSurface,
                                             freeParameters, covCache);

  variantCovariance = std::get<0>(covJacAtBound);
  variantJacobian = std::get<1>(covJacAtBound);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);
}

BOOST_AUTO_TEST_CASE(covariance_transport_curvilinear_start) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  Vector3 position{1., 2., 3.};
  Vector3 direction = Vector3(5., 6., 7.).normalized();
  ActsScalar time = 4.;
  ActsScalar qop = 0.125;

  // Free & bound parameters
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;

  auto planeSurface = Surface::makeShared<PlaneSurface>(position, direction);

  BoundSymMatrix boundCovariance = 8. * BoundSymMatrix::Identity();

  CovarianceCache covCache(position, direction, boundCovariance);
  // Test that the transport should be applied now
  BOOST_CHECK(covCache.applyTransport == true);
  // Test that the set covariance is 5x5 and what it has been set
  BOOST_CHECK_EQUAL(std::get<BoundSymMatrix>(covCache.covariance),
                    boundCovariance);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(covCache.covariance),
                    std::bad_variant_access);

  // Test that the bound to free jacobian has a value
  BOOST_CHECK(covCache.boundToFreeJacobian.has_value());
  BOOST_CHECK(not covCache.boundToFreeJacobian.value().isApprox(
      BoundToFreeMatrix::Zero()));
  // Test that the free transport jacobian, derivative is properly set up
  BOOST_CHECK_EQUAL(covCache.freeTransportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(covCache.freeToPathDerivatives, FreeVector::Zero());
  // Check that the surface is set
  BOOST_CHECK_EQUAL(covCache.atSurface, nullptr);
  BOOST_CHECK_EQUAL(covCache.atPosition, position);

  // Transport bound to the same bound surface
  auto covJacAtBound = transportCovarianceToBound(tgContext, *planeSurface,
                                                  freeParameters, covCache);
  auto variantCovariance = std::get<0>(covJacAtBound);
  auto variantJacobian = std::get<1>(covJacAtBound);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to curvilinear on the same place
  auto covJacAtCurv = transportCovarianceToCurvilinear(direction, covCache);
  variantCovariance = std::get<0>(covJacAtCurv);
  variantJacobian = std::get<1>(covJacAtCurv);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to Free on the same place
  auto covJacAtFree = transportCovarianceToFree(covCache);
  variantCovariance = std::get<0>(covJacAtFree);
  variantJacobian = std::get<1>(covJacAtFree);
  BOOST_CHECK_THROW(std::get<BoundSymMatrix>(variantCovariance),
                    std::bad_variant_access);
}

BOOST_AUTO_TEST_CASE(covariance_transport_free_start) {
  // Create a test context
  GeometryContext tgContext = GeometryContext();

  // Some start parameters
  Vector3 position{1., 2., 3.};
  ActsScalar time = 4.;
  Vector3 direction = Vector3(5., 6., 7.).normalized();

  auto planeSurface = Surface::makeShared<PlaneSurface>(position, direction);

  ActsScalar qop = 0.125;
  FreeVector freeParameters;
  freeParameters << position[0], position[1], position[2], time, direction[0],
      direction[1], direction[2], qop;

  FreeSymMatrix freeCovariance = 8. * FreeSymMatrix::Identity();

  CovarianceCache covCache(freeParameters, freeCovariance);
  BOOST_CHECK(covCache.applyTransport == true);
  // Test that the set covariance is 5x5 and what it has been set
  BOOST_CHECK_THROW(std::get<BoundSymMatrix>(covCache.covariance),
                    std::bad_variant_access);
  BOOST_CHECK_EQUAL(std::get<FreeSymMatrix>(covCache.covariance),
                    freeCovariance);
  // Test that the bound to free jacobian has NO value
  BOOST_CHECK(not covCache.boundToFreeJacobian.has_value());
  // Test that the free transport jacobian, derivative is properly set up
  BOOST_CHECK_EQUAL(covCache.freeTransportJacobian, FreeMatrix::Identity());
  BOOST_CHECK_EQUAL(covCache.freeToPathDerivatives, FreeVector::Zero());
  // Check that the surface is NOT set
  BOOST_CHECK_EQUAL(covCache.atSurface, nullptr);
  BOOST_CHECK_EQUAL(covCache.atPosition, position);

  // Transport bound to the same bound surface
  auto covJacAtBound = transportCovarianceToBound(tgContext, *planeSurface,
                                                  freeParameters, covCache);

  auto variantCovariance = std::get<0>(covJacAtBound);
  auto variantJacobian = std::get<1>(covJacAtBound);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to curvilinear on the same place
  auto covJacAtCurv = transportCovarianceToCurvilinear(direction, covCache);
  variantCovariance = std::get<0>(covJacAtCurv);
  variantJacobian = std::get<1>(covJacAtCurv);
  BOOST_CHECK_THROW(std::get<FreeSymMatrix>(variantCovariance),
                    std::bad_variant_access);

  // Transport bound to Free on the same place
  auto covJacAtFree = transportCovarianceToFree(covCache);
  variantCovariance = std::get<0>(covJacAtFree);
  variantJacobian = std::get<1>(covJacAtFree);
  BOOST_CHECK_THROW(std::get<BoundSymMatrix>(variantCovariance),
                    std::bad_variant_access);
}

}  // namespace Test
}  // namespace Acts