// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/LineWithPartials.hpp"

#include <random>

using namespace Acts;
using namespace Acts::detail;
using Line_t = LineWithPartials<double>;
using ParamVector = Line_t::ParamVector;
using ParIndices = Line_t::ParIndices;
/// The random number generator used in the framework.
using RandomEngine = std::mt19937;  ///< Mersenne Twister

ParamVector makeTrack(RandomEngine& rndEngine) {
  using namespace Acts::UnitLiterals;
  ParamVector pars{ParamVector::Zero()};
  pars[ParIndices::x0] =
      static_cast<double>(rndEngine() % 1000) - 500;  // Random x0 in [-50, 50]
  pars[ParIndices::y0] =
      static_cast<double>(rndEngine() % 1000) - 500;  // Random y0 in [-50, 50]
  pars[ParIndices::theta] = static_cast<double>(rndEngine() % 180) *
                            1_degree;  // Random theta in [0, 180)
  pars[ParIndices::phi] = static_cast<double>(rndEngine() % 360) *
                          1_degree;  // Random phi in [-180, 180)
  return pars;
}

BOOST_AUTO_TEST_SUITE(LineWithPartialTests)

BOOST_AUTO_TEST_CASE(lineParameterTest) {
  Line_t newLine{};
  constexpr unsigned trials = 1000;
  constexpr double tolerance = 1.e-12;
  RandomEngine rndEngine{3585};
  for (unsigned i = 0; i < trials; ++i) {
    auto pars = makeTrack(rndEngine);
    newLine.updateParameters(pars);
    BOOST_CHECK_CLOSE(newLine.position()[Acts::eX], pars[ParIndices::x0],
                      tolerance);
    BOOST_CHECK_CLOSE(newLine.position()[Acts::eY], pars[ParIndices::y0],
                      tolerance);
    const Acts::Vector3 dir = Acts::makeDirectionFromPhiTheta(
        pars[ParIndices::phi], pars[ParIndices::theta]);
    BOOST_CHECK_CLOSE(newLine.direction().dot(dir), 1., tolerance);
  }
}
BOOST_AUTO_TEST_CASE(lineGradientTest) {
  Line_t newLine{};
  constexpr unsigned trials = 10;
  RandomEngine rndEngine{26934};
  constexpr double h = 1.e-8;
  constexpr double tolerance = 1.e-7;
  for (unsigned trial = 0; trial < trials; ++trial) {
    const ParamVector pars{makeTrack(rndEngine)};
    std::cout << "lineGradientTest -- Generated parameters: "
              << pars.transpose() << std::endl;
    Line_t segLine{};
    segLine.updateParameters(pars);
    for (const int param : {ParIndices::theta, ParIndices::phi}) {
      ParamVector parsUp{pars}, parsDn{pars};
      parsUp[param] += h;
      parsDn[param] -= h;
      Line_t segLineUp{}, segLineDn{};
      segLineUp.updateParameters(parsUp);
      segLineDn.updateParameters(parsDn);

      const Vector3 numDeriv{(segLineUp.direction() - segLineDn.direction()) /
                             (2. * h)};
      BOOST_CHECK_LE((numDeriv - segLine.gradient(param)).norm(), tolerance);
      /** Calculate the second order derivatives of the line partials */
      for (const int param1 : {ParIndices::theta, ParIndices::phi}) {
        ParamVector parsUp1{pars}, parsDn1{pars};
        parsUp1[param1] += h;
        parsDn1[param1] -= h;

        segLineUp.updateParameters(parsUp1);
        segLineDn.updateParameters(parsDn1);

        const Vector3 numDeriv1{
            (segLineUp.gradient(param) - segLineDn.gradient(param)) / (2. * h)};
        const Vector3& analyticDeriv = segLine.hessian(param, param1);
        BOOST_CHECK_LE((numDeriv1 - analyticDeriv).norm(), tolerance);
      }
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
