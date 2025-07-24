// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/Line3DWithPartialDerivatives.hpp"

#include <random>

using namespace Acts;
using namespace Acts::detail;
using Line_t = Line3DWithPartialDerivatives<double>;
using ParamVector = Line_t::ParamVector;
using ParIndices = Line_t::ParIndices;
using namespace Acts::UnitLiterals;
/// The random number generator used in the framework.
using RandomEngine = std::mt19937;  ///< Mersenne Twister

ParamVector makeTrack(RandomEngine& rndEngine) {
  ParamVector pars{};
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
  constexpr std::size_t trials = 1000;
  constexpr double tolerance = 1.e-12;
  RandomEngine rndEngine{3585};
  for (std::size_t i = 0; i < trials; ++i) {
    auto pars = makeTrack(rndEngine);
    newLine.updateParameters(pars);
    BOOST_CHECK_CLOSE(newLine.position()[Acts::eX], pars[ParIndices::x0],
                      tolerance);
    BOOST_CHECK_CLOSE(newLine.position()[Acts::eY], pars[ParIndices::y0],
                      tolerance);
    BOOST_CHECK_LE(newLine.position()[Acts::eZ], tolerance);
    const Acts::Vector3 dir = Acts::makeDirectionFromPhiTheta(
        pars[ParIndices::phi], pars[ParIndices::theta]);
    BOOST_CHECK_CLOSE(newLine.direction().dot(dir), 1., tolerance);
  }
}
BOOST_AUTO_TEST_CASE(lineGradientTest) {
  constexpr std::size_t trials = 1000;
  RandomEngine rndEngine{26934};
  constexpr double h = 1.e-8;
  constexpr double tolerance = 1.e-7;
  for (std::size_t trial = 0; trial < trials; ++trial) {
    const ParamVector pars{makeTrack(rndEngine)};
    std::cout << "lineGradientTest -- Generated parameters  x: "
              << pars[ParIndices::x0] << ", y: " << pars[ParIndices::y0]
              << ", theta: " << (pars[ParIndices::theta] / 1_degree)
              << ", phi: " << (pars[ParIndices::phi] / 1_degree) << std::endl;
    Line_t segLine{};
    segLine.updateParameters(pars);

    BOOST_CHECK_LE((segLine.gradient(ParIndices::x0) - Vector3::UnitX()).norm(),
                   tolerance);
    BOOST_CHECK_LE((segLine.gradient(ParIndices::y0) - Vector3::UnitY()).norm(),
                   tolerance);

    for (const std::size_t param : {ParIndices::theta, ParIndices::phi}) {
      break;
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
      for (const std::size_t param1 : {ParIndices::theta, ParIndices::phi,
                                       ParIndices::x0, ParIndices::y0}) {
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
