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
/// @brief Express a time in terms of nanoseconds
constexpr double inNanoS(const double x) {
  return x / 1._ns;
}

}  // namespace
namespace Acts::Experimental::detail {

void FastStrawLineFitter::FitAuxiliaries::print(std::ostream& ostr) const {
  ostr << std::format("T_zzyy: {:.3f}, ", T_zzyy)
       << std::format("T_yz: {:.3f}, ", T_yz)
       << std::format("T_rz: {:.3f}, ", T_rz)
       << std::format("T_ry: {:.3f}, ", T_ry)
       << std::format("centre ( {:.3f}, {:3f}), ", centerY, centerZ)
       << "y0: " << fitY0 << ", norm: " << covNorm << ", nDoF: " << nDoF;
}
void FastStrawLineFitter::FitAuxiliariesWithT0 ::print(
    std::ostream& ostr) const {
  FitAuxiliaries::print(ostr);
  ostr << std::format("T_vz: {:.3f}, ", T_vz)
       << std::format("T_vy: {:.3f}, ", T_vy)
       << std::format("T_az: {:.3f}, ", T_az)
       << std::format("T_ay: {:.3f}, ", T_ay)
       << std::format("R_vr: {:.3f}, ", R_vr)
       << std::format("R_va: {:.3f}, ", R_va)
       << std::format("R_vv: {:.3f}, ", R_vv) << " --- "
       << std::format("y_{{0}}^{{''}}:  {:.3f}, ", fitY0TwoPrime)
       << std::format("y_{{0}}^{{'}} :  {:.3f}", fitY0Prime);
}
using Vector = FastStrawLineFitter::Vector;
const Acts::Logger& FastStrawLineFitter::logger() const {
  assert(m_logger != nullptr);
  return *m_logger;
}
FastStrawLineFitter::FastStrawLineFitter(const Config& cfg,
                                         std::unique_ptr<const Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {}

void FastStrawLineFitter::calcAngularDerivatives(const TrigonomHelper& angles,
                                                 const FitAuxiliaries& fitPars,
                                                 double& thetaPrime,
                                                 double& thetaTwoPrime) const {
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
std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const FitAuxiliaries& fitPars) const {
  /// No degrees of freedom -> no valid parameters
  if (!fitPars.nDoF) {
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Input fit parameters "
                      << fitPars);

  const double thetaGuess = startTheta(fitPars);
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
    const TrigonomHelper angles{result.theta};
    calcAngularDerivatives(angles, fitPars, thetaPrime, thetaTwoPrime);
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_VERBOSE(
        __func__ << "() - " << __LINE__ << ": Fit iteration #" << result.nIter
                 << std::format("-- theta: {:.3f}, ", inDeg(result.theta))
                 << std::format(" thetaPrime: {:.3f}, ", inDeg(thetaPrime))
                 << std::format(" thetaTwoPrime: {:.3f}", inDeg(thetaTwoPrime))
                 << " --> step-size: " << std::format("{:.8f}", inDeg(update))
                 << std::format(
                        " --> next: {:.3f} ",
                        inDeg(result.theta - thetaPrime / thetaTwoPrime)));

    if (std::abs(update) < m_cfg.precCutOff) {
      converged = true;
    }
    result.theta -= update;
  }

  if (!converged) {
    ACTS_WARNING(
        __func__ << "() - " << __LINE__
                 << ": The fast straw fit did not converge " << fitPars
                 << ", \n"
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
  ACTS_INFO(
      __func__ << "() - " << __LINE__ << ": Fit converged in #" << result.nIter
               << " iterations with final parameters: "
               << std::format("theta: {:.3f} pm {:.3f}, ", inDeg(result.theta),
                              inDeg(result.dTheta))
               << std::format("y0: {:.3f}  pm {:.3f}", result.y0, result.dY0));
  return result;
}
#ifdef STONJEK
std::optional<FastStrawLineFitter::FitResultT0> FastStrawLineFitter::fit(
    const FitAuxiliariesWithT0& fitPars) const {
  /// No degrees of freedom -> no valid parameters
  if (!fitPars.nDoF) {
    return std::nullopt;
  }
  auto noT0Result = fit(static_cast<const FitAuxiliaries&>(fitPars));

  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Estimated parameters "
                     << fitPars);

  FitResultT0 result{};
  result.nDoF = fitPars.nDoF;

  ActsSquareMatrix<2> cov{ActsSquareMatrix<2>::Zero()};
  Vector2 grad{Vector2::Zero()};
  Vector2 pars{Vector2::Zero()};

  pars[0] = startTheta(fitPars);
  pars[1] = 82._ns;
  bool converged{false};
  while (converged == false && result.nIter <= m_cfg.maxIter) {
    ++result.nIter;
    const TrigonomHelper angles{pars[0]};
    calcAngularDerivatives(angles, fitPars, grad[0], cov(0, 0));

    cov(1, 0) = cov(0, 1) =
        fitPars.T_vz * angles.cosTheta + fitPars.T_vy * angles.sinTheta;
    cov(1, 1) = -fitPars.fitY0Prime * fitPars.fitY0Prime -
                fitPars.fitY0Prime * fitPars.fitY0TwoPrime -
                fitPars.T_az * angles.sinTheta -
                fitPars.T_ay * angles.cosTheta + fitPars.R_vv + fitPars.R_va;

    grad[1] = fitPars.fitY0 * fitPars.fitY0Prime - fitPars.R_vr +
              fitPars.T_vz * angles.sinTheta - fitPars.T_vy * angles.cosTheta;

    const Vector2 update = cov.inverse() * grad;
    if (update.norm() < m_cfg.precCutOff) {
      converged = true;
      break;
    }
    ACTS_INFO(__func__ << "() - " << __LINE__ << " Iteration: " << result.nIter
                       << ", theta: " << inDeg(pars[0])
                       << ", time: " << inNanoS(pars[1]) << " gradient: ("
                       << inDeg(grad[0]) << ", " << inNanoS(grad[1])
                       << "), covariance:" << std::endl
                       << toString(cov) << std::endl
                       << " update: (" << inDeg(update[0]) << ", "
                       << inNanoS(update[1]) << ").");
    pars -= update;
  }
  return std::nullopt;
}
#endif
}  // namespace Acts::Experimental::detail
