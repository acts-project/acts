// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/CovarianceTransport.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <random>

int main(int argc, char* argv[]) {
  using namespace Acts;

  size_t iterations = 3;
  size_t runs = 1000;
  if (argc >= 2) {
    iterations = std::stoi(argv[1]);
  }
  if (argc >= 3) {
    runs = std::stoi(argv[2]);
  }

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

  std::minstd_rand rng;
  std::uniform_real_distribution<> uniform(0.5, 0.95);

  unsigned int sillyCounter = 0;

  ACTS_LOCAL_LOGGER(
      getDefaultLogger("CovarianceTransport", Acts::Logging::Level(0)));

  const auto cov_transport_bound_bound = Acts::Test::microBenchmark(
      [&] {
        BoundSquareMatrix boundCovariance =
            uniform(rng) * BoundSquareMatrix::Identity();
        boundCovariance(eBoundLoc0, eBoundPhi) = 0.076;
        boundCovariance(eBoundPhi, eBoundLoc0) = 0.076;
        boundCovariance(eBoundLoc0, eBoundQOverP) = -0.022;
        boundCovariance(eBoundQOverP, eBoundLoc0) = -0.022;
        boundCovariance(eBoundLoc1, eBoundTheta) = -0.007;
        boundCovariance(eBoundTheta, eBoundLoc1) = -0.007;

        CovarianceCache covCache(tgContext, *planeSurface, position,
                                 boundParameters, boundCovariance);

        // Transport bound to another bound surface
        auto covJacAtBound = transportCovarianceToBound(
            tgContext, *otherSurface, freeParameters, covCache);

        // Mimic access to the result
        const auto& variantCovariance = std::get<0>(covJacAtBound);
        const auto& variantJacobian = std::get<1>(covJacAtBound);

        const auto& covariance = std::get<BoundSquareMatrix>(variantCovariance);
        const auto& jacobian = std::get<BoundMatrix>(variantJacobian);

        if (covariance(eBoundLoc0, eBoundLoc0) > 0 and
            jacobian(eBoundLoc0, eBoundLoc1) > 0.) {
          ++sillyCounter;
        }
      },
      iterations, runs);

  ACTS_INFO("Execution stats: " << cov_transport_bound_bound);
}
