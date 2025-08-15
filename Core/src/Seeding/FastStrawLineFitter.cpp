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
  /// No degrees of freedom -> no valid parameters
  if (!fitPars.nDoF) {
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Input fit parameters "
                      << std::format("T_zzyy: {:.3f}, ", fitPars.T_zzyy)
                      << std::format("T_yz: {:.3f}, ", fitPars.T_yz)
                      << std::format("T_rz: {:.3f}, ", fitPars.T_rz)
                      << std::format("T_ry: {:.3f}, ", fitPars.T_ry)
                      << std::format("centre ( {:.3f}, {:3f}), ",
                                     fitPars.centerY, fitPars.centerZ)
                      << "y0: " << fitPars.fitY0 << ", norm: "
                      << fitPars.covNorm << ", nDoF: " << fitPars.nDoF);

  double thetaGuess =
      std::atan2(2. * (fitPars.T_yz - fitPars.T_ry), fitPars.T_zzyy) / 2.;
  if (thetaGuess < 0) {
    thetaGuess += std::numbers::pi;
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__
                      << ": Start fast fit seed  theta: " << inDeg(thetaGuess));
  ////
  bool converged{false};
  FitResult result{};
  result.theta = thetaGuess;

  result.nDoF = fitPars.nDoF;
  double thetaPrime{0.}, thetaTwoPrime{0.};
  while ((converged == false) && result.nIter <= m_cfg.maxIter) {
    ++result.nIter;
    const double cosTheta = std::cos(result.theta);
    const double sinTheta = std::sin(result.theta);
    const double cosTwoTheta = std::cos(2. * result.theta);
    const double sinTwoTheta = std::sin(2. * result.theta);
    thetaPrime = 0.5 * fitPars.T_zzyy * sinTwoTheta -
                 fitPars.T_yz * cosTwoTheta - fitPars.T_rz * cosTheta -
                 fitPars.T_ry * sinTheta;
    thetaTwoPrime = fitPars.T_zzyy * cosTwoTheta +
                    2. * fitPars.T_yz * sinTwoTheta + fitPars.T_rz * sinTheta -
                    fitPars.T_ry * cosTheta;
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_VERBOSE(
        __func__ << "() - " << __LINE__ << ": Fit iteration #" << result.nIter
                 << std::format("-- theta: {:.3f}, ", inDeg(result.theta))
                 << std::format(" thetaPrime: {:.3f}, ", inDeg(thetaPrime))
                 << std::format(" thetaTwoPrime: {:.3f}", inDeg(thetaTwoPrime))
                 << " --> step-size: " << std::format("{:.8f}", inDeg(update))
                 << std::format(
                        " --> next theta: {:.3f} ",
                        inDeg(result.theta - thetaPrime / thetaTwoPrime)));

    if (std::abs(update) < m_cfg.precCutOff) {
      converged = true;
    }
    result.theta -= update;
  }

  if (!converged) {
    ACTS_WARNING(
        __func__ << "() - " << __LINE__
                 << ": The fast straw fit did not converge nDoF: "
                 << fitPars.nDoF << ", "
                 << std::format("-- theta: {:.3f}, ", inDeg(result.theta))
                 << std::format(" thetaPrime: {:.3f}, ", inDeg(thetaPrime))
                 << std::format(" thetaTwoPrime: {:.3f}", inDeg(thetaTwoPrime))
                 << " --> "
                 << std::format("{:.8f}", inDeg(thetaPrime / thetaTwoPrime)));
    return std::nullopt;
  }
  result.dTheta = std::sqrt(1. / thetaTwoPrime);
  const double tanTheta = std::tan(result.theta);
  const double secTheta = 1. / std::cos(result.theta);
  result.y0 =
      fitPars.centerY - fitPars.centerZ * tanTheta + fitPars.fitY0 * secTheta;
  result.dY0 = Acts::fastHypot(-fitPars.centerZ * Acts::pow(secTheta, 2) +
                                   fitPars.fitY0 * secTheta * tanTheta,
                               secTheta);

  result.dY0 *= result.dTheta;
  ACTS_DEBUG(
      __func__ << "() - " << __LINE__ << ": Fit converged in #" << result.nIter
               << " iterations with final parameters: "
               << std::format("theta: {:.3f} pm {:.3f}, ", inDeg(result.theta),
                              inDeg(result.dTheta))
               << std::format("y0: {:.3f}  pm {:.3f}", result.y0, result.dY0));
  return result;
}

std::optional<FastStrawLineFitter::FitResultT0> FastStrawLineFitter::fit(
    const FitAuxiliariesWithT0& fitPars) const {
  /// No degrees of freedom -> no valid parameters
  if (!fitPars.nDoF) {
    return std::nullopt;
  }

  ACTS_INFO(
      __func__ << "() - " << __LINE__ << ": Estimated "
               << std::format("T_zzyy: {:.3f}, ", fitPars.T_zzyy)
               << std::format("T_yz: {:.3f}, ", fitPars.T_yz)
               << std::format("T_rz: {:.3f}, ", fitPars.T_rz)
               << std::format("T_ry: {:.3f}, ", fitPars.T_ry)
               << std::format("centre ( {:.3f}, {:3f}), ", fitPars.centerY,
                              fitPars.centerZ)
               << "y0: " << fitPars.fitY0 << ", norm: " << fitPars.covNorm
               << ", nDoF: " << fitPars.nDoF << "\n"
               << std::format("T_vz: {:.3f}, ", fitPars.T_vz)
               << std::format("T_vy: {:.3f}, ", fitPars.T_vy)
               << std::format("T_az: {:.3f}, ", fitPars.T_az)
               << std::format("T_ay: {:.3f}, ", fitPars.T_ay)
               << std::format("R_vr: {:.3f}, ", fitPars.R_vr)
               << std::format("R_va: {:.3f}, ", fitPars.R_va)
               << std::format("R_vv: {:.3f}, ", fitPars.R_vv) << " --- "
               << std::format("y_{{0}}^{{''}}:  {:.3f}, ",
                              fitPars.fitY0TwoPrime)
               << std::format("y_{{0}}^{{'}} :  {:.3f}", fitPars.fitY0Prime));

  FitResultT0 result{};
  return std::nullopt;
}

}  // namespace Acts::Experimental::detail
