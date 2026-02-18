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
#include "Acts/Utilities/detail/periodic.hpp"

#include <format>

namespace {
using namespace Acts::UnitLiterals;
/// @brief Expresses an angle in degree
/// @param x: Angle in radians
constexpr double inDeg(const double x) {
  return x / 1._degree;
}
/// @brief Express a time in units of nanoseconds
constexpr double inNanoS(const double x) {
  return x / 1._ns;
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
void FastStrawLineFitter::FitAuxiliariesWithT0 ::print(
    std::ostream& ostr) const {
  ostr << std::format("R_vr: {:.3f}, ", R_vr);
  ostr << std::format("R_ar: {:.3f}, ", R_ar);
  ostr << std::format("R_vv: {:.3f}, ", R_vv);
  ostr << std::format("R_a: {:.3f}, ", R_a);
  ostr << std::format("R_v: {:.3f}, ", R_v);
  ostr << std::format("T_vz: {:.3f}, ", T_vz);
  ostr << std::format("T_vy: {:.3f}, ", T_vy);
  ostr << std::format("T_az: {:.3f}, ", T_az);
  ostr << std::format("T_ay: {:.3f}, ", T_ay);
  FitAuxiliaries::print(ostr);
}

void FastStrawLineFitter::FitResult::print(std::ostream& ostr) const {
  ostr << "# iteration: " << nIter << ", nDoF: " << nDoF << ", chi2: " << chi2
       << ", chi2 / nDoF: " << (chi2 / static_cast<double>(nDoF)) << ",\n";
  ostr << std::format("theta: {:.3f} pm {:.3f}, ", inDeg(theta), inDeg(dTheta))
       << std::format("y0: {:.3f}  pm {:.3f}", y0, dY0);
}
void FastStrawLineFitter::FitResultT0::print(std::ostream& ostr) const {
  FitResult::print(ostr);
  ostr << std::format(", t0: {:.3f} pm {:.3f}", inNanoS(t0), inNanoS(dT0));
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
  return Acts::detail::wrap_periodic(thetaGuess, 0., std::numbers::pi);
}
double FastStrawLineFitter::calcTimeGrad(const TrigonomHelper& angles,
                                         const FitAuxiliariesWithT0& pars) {
  return pars.fitY0 * pars.R_v - pars.R_vr + pars.T_vz * angles.sinTheta -
         pars.T_vy * angles.cosTheta;
}

void FastStrawLineFitter::completeResult(const FitAuxiliaries& fitPars,
                                         FitResult& result) const {
  const double tanTheta = std::tan(result.theta);
  const double secTheta = 1. / std::cos(result.theta);
  result.y0 =
      fitPars.centerY - fitPars.centerZ * tanTheta + fitPars.fitY0 * secTheta;
  result.dY0 = Acts::fastHypot(-fitPars.centerZ * Acts::square(secTheta) +
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
    if (std::abs(thetaTwoPrime) < 2. * std::numeric_limits<double>::epsilon()) {
      ACTS_WARNING(__func__
                   << "() - " << __LINE__
                   << ": The fast straw fit encountered a singular "
                      "second derivative "
                   << fitPars << ", \n"
                   << printThetaStep(result.theta, thetaPrime, thetaTwoPrime)
                   << "\n"
                   << result);
      return std::nullopt;
    }
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_VERBOSE(
        __func__ << "() - " << __LINE__ << ": Fit iteration #" << result.nIter
                 << " "
                 << printThetaStep(result.theta, thetaPrime, thetaTwoPrime));

    if (std::abs(update) < m_cfg.precCutOff) {
      result.dTheta = std::sqrt(1. / Acts::abs(thetaTwoPrime));
      completeResult(fitPars, result);
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

FastStrawLineFitter::UpdateStatus FastStrawLineFitter::updateIteration(
    const FitAuxiliariesWithT0& fitPars, FitResultT0& fitResult) const {
  ++fitResult.nIter;
  if (fitResult.nIter > m_cfg.maxIter) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__
                        << ": The fast straw fit did not converge " << fitPars
                        << "\n"
                        << fitResult);
    return UpdateStatus::Exceeded;
  }

  UpdateStatus retCode{UpdateStatus::GoodStep};
  SquareMatrix<2> cov{SquareMatrix<2>::Zero()};
  Vector2 grad{Vector2::Zero()};

  const TrigonomHelper angles{fitResult.theta};

  calcAngularDerivatives(angles, fitPars, grad[0], cov(0, 0));
  grad[1] = calcTimeGrad(angles, fitPars);
  if (grad.norm() < m_cfg.precCutOff) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Fit converged "
                        << fitResult);
    retCode = UpdateStatus::Converged;
  }

  cov(1, 0) = cov(0, 1) =
      fitPars.T_vz * angles.cosTheta + fitPars.T_vy * angles.sinTheta;

  cov(1, 1) = -fitPars.R_v * fitPars.R_v * fitPars.covNorm -
              fitPars.fitY0 * fitPars.R_a + fitPars.R_vv + fitPars.R_ar -
              (fitPars.T_az * angles.sinTheta - fitPars.T_ay * angles.cosTheta);
  ACTS_VERBOSE(
      __func__
      << "() - " << __LINE__
      << ": -fitPars.R_v * fitPars.R_v * fitPars.covNorm + fitPars.R_vv: "
      << (-fitPars.R_v * fitPars.R_v * fitPars.covNorm + fitPars.R_vv)
      << ", fitPars.fitY0 * fitPars.R_a: " << (fitPars.fitY0 * fitPars.R_a)
      << ", fitPars.R_ar: " << fitPars.R_ar
      << ", fitPars.T_az * angles.sinTheta - fitPars.T_ay * angles.cosTheta: "
      << (fitPars.T_az * angles.sinTheta - fitPars.T_ay * angles.cosTheta));
  const auto invCov = cov.inverse();
  if (invCov(1, 1) < 0) {
    ACTS_DEBUG("Invalid covariance\n"
               << invCov << cov << ", determinant: " << cov.determinant()
               << ", " << fitPars);
    return UpdateStatus::Exceeded;
  }
  const Vector2 update = invCov * grad;
  // We compute also the normalized update, defined as the parameter
  // update expressed in units of the parameter uncertainties. This quantifies
  // the significance of the update relative to the estimated errors.
  double normUpdate{0.};
  for (unsigned i = 0; i < 2; ++i) {
    normUpdate += Acts::square(update[i]) / invCov(i, i);
  }

  ACTS_VERBOSE(__func__ << "() - " << __LINE__ << " intermediate result "
                        << fitResult << "\n"
                        << std::format("gradient: ({:.3f}, {:.3f})",
                                       inDeg(grad[0]), inNanoS(grad[1]))
                        << ", covariance:" << std::endl
                        << toString(cov) << std::endl
                        << cov.determinant()
                        << std::format(" update: ({:.3f}, {:.3f}),",
                                       inDeg(update[0]), inNanoS(update[1]))
                        << " normUpdate: " << std::sqrt(normUpdate));

  if (std::sqrt(normUpdate) < m_cfg.normPrecCutOff ||
      update.norm() < m_cfg.precCutOff) {
    ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Fit converged "
                        << fitResult);
    retCode = UpdateStatus::Converged;
  }

  if (retCode == UpdateStatus::Converged) {
    fitResult.dTheta = std::sqrt(invCov(0, 0));
    fitResult.dT0 = std::sqrt(invCov(1, 1));
    completeResult(fitPars, fitResult);
    return retCode;
  }

  fitResult.t0 -= update[1];
  fitResult.theta -= update[0];
  return retCode;
}
}  // namespace Acts::Experimental::detail
