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
/// @brief Expresses an angle in degree
/// @param x: Angle in radians
constexpr double inDeg(const double x) {
  return x / 1._degree;
}
std::string printThetaStep(const double theta, const double thetaPrime,
                           const double thetaTwoPrime) {
  return std::format(
      "-- theta: {:.3f}, thetaPrime: {:.3f}, thetaTwoPrime: {:.3f}  --> "
      "step-size: {:.8f}",
      inDeg(theta), inDeg(thetaPrime), inDeg(thetaTwoPrime),
      inDeg(thetaPrime / thetaTwoPrime));
}
}  // namespace
namespace Acts::Experimental::detail {

void FastStrawLineFitter::FitAuxiliaries::print(std::ostream& ostr) const {
  ostr << std::format("T_zzyy: {:.3f}, ", T_zzyy)
       << std::format("T_yz: {:.3f}, ", T_yz)
       << std::format("T_rz: {:.3f}, ", T_rz)
       << std::format("T_ry: {:.3f}, ", T_ry)
       << std::format("centre ( {:.3f}, {:3f}), ", centerY, centerZ)
       << std::format("y0: {:.3f}", fitY0)
       << std::format(", norm: {:.3f}", covNorm) << ", nDoF: " << nDoF;
}
void FastStrawLineFitter::FitResult::print(std::ostream& ostr) const {
  ostr << "# iteration: " << nIter << ", nDoF: " << nDoF << ", chi2: " << chi2
       << ", chi2 / nDoF: " << (chi2 / static_cast<double>(nDoF)) << ",\n";
  ostr << std::format("theta: {:.3f} pm {:.3f}, ", inDeg(theta), inDeg(dTheta))
       << std::format("y0: {:.3f}  pm {:.3f}", y0, dY0);
}

using Vector = FastStrawLineFitter::Vector;
const Acts::Logger& FastStrawLineFitter::logger() const {
  assert(m_logger != nullptr);
  return *m_logger;
}
FastStrawLineFitter::FastStrawLineFitter(const Config& cfg,
                                         std::unique_ptr<const Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {}

inline void FastStrawLineFitter::calcAngularDerivatives(
    const TrigonomHelper& angles, const FitAuxiliaries& fitPars,
    double& thetaPrime, double& thetaTwoPrime) const {
  thetaPrime = 0.5 * fitPars.T_zzyy * angles.sinTwoTheta -
               fitPars.T_yz * angles.cosTwoTheta -
               fitPars.T_rz * angles.cosTheta - fitPars.T_ry * angles.sinTheta;
  thetaTwoPrime = fitPars.T_zzyy * angles.cosTwoTheta +
                  2. * fitPars.T_yz * angles.sinTwoTheta +
                  fitPars.T_rz * angles.sinTheta -
                  fitPars.T_ry * angles.cosTheta;
}
double FastStrawLineFitter::startTheta(const FitAuxiliaries& fitPars) {
  double thetaGuess =
      std::atan2(2. * (fitPars.T_yz - fitPars.T_ry), fitPars.T_zzyy) / 2.;
  return thetaGuess + (thetaGuess < 0 ? std::numbers::pi : 0.);
}
void FastStrawLineFitter::completeResult(const FitAuxiliaries& fitPars,
                                         const double thetaTwoPrime,
                                         FitResult& result) const {
  result.dTheta = std::sqrt(1. / thetaTwoPrime);
  const double tanTheta = std::tan(result.theta);
  const double secTheta = 1. / std::cos(result.theta);
  result.y0 =
      fitPars.centerY - fitPars.centerZ * tanTheta + fitPars.fitY0 * secTheta;
  result.dY0 = Acts::fastHypot(-fitPars.centerZ * Acts::pow(secTheta, 2) +
                                   fitPars.fitY0 * secTheta * tanTheta,
                               secTheta);

  result.dY0 *= result.dTheta;
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Fit succeeded " << result);
}

std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const FitAuxiliaries& fitPars) const {
  /// No degrees of freedom -> no valid parameters
  if (fitPars.nDoF == 0u) {
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Input fit parameters "
                      << fitPars);

  const double thetaGuess = startTheta(fitPars);
  ACTS_DEBUG(__func__ << "() - " << __LINE__
                      << ": Start fast fit seed  theta: " << inDeg(thetaGuess));
  ////
  FitResult result{};
  result.theta = thetaGuess;
  result.nDoF = fitPars.nDoF;

  double thetaPrime{0.};
  double thetaTwoPrime{0.};
  while (result.nIter <= m_cfg.maxIter) {
    ++result.nIter;
    const TrigonomHelper angles{result.theta};
    calcAngularDerivatives(angles, fitPars, thetaPrime, thetaTwoPrime);
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_VERBOSE(
        __func__ << "() - " << __LINE__ << ": Fit iteration #" << result.nIter
                 << " "
                 << printThetaStep(result.theta, thetaPrime, thetaTwoPrime));

    if (std::abs(update) < m_cfg.precCutOff) {
      completeResult(fitPars, thetaTwoPrime, result);
      return result;
    }
    result.theta -= update;
  }

  ACTS_WARNING(
      __func__ << "() - " << __LINE__
               << ": The fast straw fit did not converge " << fitPars << ", \n"
               << printThetaStep(result.theta, thetaPrime, thetaTwoPrime)
               << "\n"
               << result);
  return std::nullopt;
}

}  // namespace Acts::Experimental::detail
