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
       << std::format("y0: {:.3f}", fitY0)
       << std::format(", norm: {:.3f}", covNorm) << ", nDoF: " << nDoF;
}
void FastStrawLineFitter::FitAuxiliariesWithT0 ::print(
    std::ostream& ostr) const {
  FitAuxiliaries::print(ostr);
  ostr << ", " << std::format("T_vz: {:.3f}, ", T_vz)
       << std::format("T_vy: {:.3f}, ", T_vy)
       << std::format("T_az: {:.3f}, ", T_az)
       << std::format("T_ay: {:.3f}, ", T_ay)
       << std::format("R_vr: {:.3f}, ", R_vr)
       << std::format("R_va: {:.3f}, ", R_va)
       << std::format("R_vv: {:.3f}, ", R_vv) << " --- "
       << std::format("y_{{0}}^{{''}}:  {:.3f}, ", fitY0TwoPrime)
       << std::format("y_{{0}}^{{'}} :  {:.3f}", fitY0Prime);
}

void FastStrawLineFitter::FitResult::print(std::ostream& ostr) const {
  ostr << "# iteration: " << nIter << ", nDoF: " << nDoF << ", chi2: " << chi2
       << ", chi2 / nDoF: " << (chi2 / nDoF) << ",\n";
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
  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Input fit parameters "
                     << fitPars);

  const double thetaGuess = startTheta(fitPars);
  ACTS_DEBUG(__func__ << "() - " << __LINE__
                      << ": Start fast fit seed  theta: " << inDeg(thetaGuess));
  ////
  bool converged{false};
  FitResult result{};
  result.theta = thetaGuess;
  result.theta = 45._degree;

  result.nDoF = fitPars.nDoF;
  double thetaPrime{0.}, thetaTwoPrime{0.};

  while ((converged == false) && result.nIter <= m_cfg.maxIter) {
    ++result.nIter;
    const TrigonomHelper angles{result.theta};
    calcAngularDerivatives(angles, fitPars, thetaPrime, thetaTwoPrime);
    const double update = thetaPrime / thetaTwoPrime;
    ACTS_INFO(
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
                 << "\n"
                 << result);
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
  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Fit succeeded " << result);
  return result;
}
#ifdef STONJEK
std::optional<FastStrawLineFitter::FitResultT0> FastStrawLineFitter::fit(
    const FitAuxiliariesWithT0& fitPars) const {
  /// No degrees of freedom -> no valid parameters

  auto noT0Result = fit(static_cast<const FitAuxiliaries&>(fitPars));

  ACTS_INFO(__func__ << "() - " << __LINE__ << ": Estimated parameters "
                     << fitPars);

  return std::nullopt;
}
#endif
}  // namespace Acts::Experimental::detail
