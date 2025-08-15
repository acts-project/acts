// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>

using namespace Acts::UnitLiterals;

namespace {
constexpr double inDeg(const double x) {
  return x / 1._degree;
}
}  // namespace
namespace Acts::Experimental::detail {

using Vector = FastStrawLineFitter::Vector;
const Acts::Logger& FastStrawLineFitter::logger() const {
  assert(m_logger != nullptr);
  return *m_logger;
}
FastStrawLineFitter::FastStrawLineFitter(const Config& cfg,
                                         std::unique_ptr<const Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {}

std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const FitAuxiliaries& fitPars) const {
  ACTS_INFO("Estimated T_zzyy: "
            << fitPars.T_zzyy << ", T_yz: " << fitPars.T_yz
            << ", T_rz: " << fitPars.T_rz << ", T_ry: " << fitPars.T_ry
            << ", centre " << toString(fitPars.centerOfGrav)
            << ", y0: " << fitPars.fitY0 << ", norm: " << fitPars.covNorm);

  const double thetaGuess =
      std::atan2(2. * (fitPars.T_yz - fitPars.T_rz), fitPars.T_zzyy) / 2.;

  ACTS_INFO("Start fast fit seed  theta: " << inDeg(thetaGuess));
  ////
  bool converged{false};
  FitResult result{};
  result.theta = thetaGuess;
  double thetaPrime{0.}, thetaTwoPrime{0.};
  while (!converged && result.nIter++ <= m_cfg.maxIter) {
    const double cosTheta = std::cos(result.theta);
    const double sinTheta = std::sin(result.theta);
    const double cosTwoTheta = std::cos(2. * result.theta);
    const double sinTwoTheta = std::sin(2. * result.theta);
    thetaPrime = 0.5 * fitPars.T_zzyy * sinTwoTheta -
                 fitPars.T_yz * cosTwoTheta - fitPars.T_rz * cosTheta -
                 fitPars.T_ry * sinTheta;
    if (std::abs(thetaPrime) < m_cfg.precCutOff) {
      converged = true;
      break;
    }
    thetaTwoPrime = fitPars.T_zzyy * cosTwoTheta +
                    2. * fitPars.T_yz * sinTwoTheta + fitPars.T_rz * sinTheta -
                    fitPars.T_ry * cosTheta;
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_INFO("Fit iteration #"
              << result.nIter << " -- theta: " << inDeg(result.theta)
              << ", thetaPrime: " << inDeg(thetaPrime)
              << ", thetaTwoPrime: " << inDeg(thetaTwoPrime) << " -- "
              << std::format("{:.8f}", inDeg(update)) << " --> next theta "
              << inDeg(result.theta - thetaPrime / thetaTwoPrime));

    if (std::abs(update) < m_cfg.precCutOff) {
      converged = true;
      break;
    }
    result.theta -= update;
  }

  if (!converged) {
    ACTS_INFO("The fit did not converge");
    return std::nullopt;
  }
  result.dTheta = std::sqrt(1. / thetaTwoPrime);
  result.y0 = fitPars.centerOfGrav.y() -
              fitPars.centerOfGrav.z() * std::tan(result.theta) +
              fitPars.fitY0 / std::cos(result.theta);

  ACTS_INFO("Fit converged in #"
            << result.nIter << " iterations with final parameters: "
            << std::format("theta: {:.3f} pm {:.3f}, ", inDeg(result.theta),
                           inDeg(result.dTheta))
            << std::format("y0: {:.3f}  pm {:.3f}", result.y0, result.dY0));
  return result;
}

}  // namespace Acts::Experimental::detail
