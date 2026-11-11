// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace Acts::Experimental {

CompositeSpacePointLineFitter::CompositeSpacePointLineFitter(
    const Config& cfg, std::unique_ptr<const Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {
  std::ranges::sort(m_cfg.parsToUse);
}

detail::FastStrawLineFitter::Config
CompositeSpacePointLineFitter::fastFitterCfg() const {
  detail::FastStrawLineFitter::Config fastCfg{};
  fastCfg.maxIter = m_cfg.maxIter;
  fastCfg.precCutOff = m_cfg.precCutOff;
  return fastCfg;
}

void CompositeSpacePointLineFitter::toChamberFrame(
    FastFitter_t::FitResult& result, const Vector_t& eY, const Vector_t& eZ,
    const CoordinateIndices coord) {
  assert(coord == CoordinateIndices::eY || coord == CoordinateIndices::eX);
  const auto iCoord = toUnderlying(coord);
  const double sinA = std::sin(result.theta);
  const double cosA = std::cos(result.theta);
  const auto dir = copySign<Vector_t, double>(sinA * eY + cosA * eZ, sinA);
  assert(std::abs(dir.z()) > std::numeric_limits<double>::epsilon());
  const double slope = dir[iCoord] / dir.z();
  auto isect =
      PlanarHelper::intersectPlane(result.y0 * eY, dir, Vector_t::UnitZ(), 0.);
  assert(isect.isValid());
  // Derivative of the direction w.r.t. the angle, with the same sign as dir
  const auto dDir = copySign<Vector_t, double>(cosA * eY - sinA * eZ, sinA);
  const double dSlope =
      (dDir[iCoord] * dir.z() - dir[iCoord] * dDir.z()) / Acts::square(dir.z());
  const double dPosdY0 = eY[iCoord] - eY.z() * slope;
  const double dPosdA = -result.y0 * eY.z() * dSlope;
  // Jacobian ((c, slope) w.r.t. (y0', angle)) = ((dPosdY0, dPosdA), (0,
  // dSlope))
  const double varY0 = Acts::square(result.dY0);
  const double varA = Acts::square(result.dTheta);
  const double covY0A = result.covY0Theta;
  result.y0 = isect.position()[iCoord];
  result.theta = slope;
  result.dY0 =
      std::sqrt(Acts::square(dPosdY0) * varY0 + 2. * dPosdY0 * dPosdA * covY0A +
                Acts::square(dPosdA) * varA);
  result.covY0Theta = (dPosdY0 * covY0A + dPosdA * varA) * dSlope;
  result.dTheta = std::abs(dSlope) * result.dTheta;
}

void CompositeSpacePointLineFitter::fastFitCovariance(
    const FastFitter_t::FitResult& planeResult, const CoordinateIndices coord,
    const Vector_t& lineDir, CovMat_t& cov) {
  using enum FitParIndex;
  assert(coord == CoordinateIndices::eY || coord == CoordinateIndices::eX);
  const bool bending = coord == CoordinateIndices::eY;
  const auto iPos = toUnderlying(bending ? y0 : x0);
  const auto iAng = toUnderlying(bending ? theta : phi);

  const double dirTheta = VectorHelpers::theta(lineDir);
  const double sinPhi = std::sin(VectorHelpers::phi(lineDir));
  if (std::abs(sinPhi) < std::numeric_limits<double>::epsilon()) {
    return;
  }
  /// Jacobian of the angle w.r.t. the slope
  const double dAng = bending ? Acts::square(std::cos(dirTheta)) / sinPhi
                              : -1. / (std::tan(dirTheta) * sinPhi);
  cov(iPos, iPos) = Acts::square(planeResult.dY0);
  cov(iAng, iAng) = Acts::square(dAng * planeResult.dTheta);
  cov(iPos, iAng) = cov(iAng, iPos) = dAng * planeResult.covY0Theta;
}

void CompositeSpacePointLineFitter::fastFitCovariance(
    const FastFitter_t::FitResult& precResult,
    const FastFitter_t::FitResult& nonPrecResult, const Vector_t& lineDir,
    CovMat_t& cov) {
  using enum FitParIndex;
  const auto iY0 = toUnderlying(y0);
  const auto iTh = toUnderlying(theta);
  const auto iX0 = toUnderlying(x0);
  const auto iPh = toUnderlying(phi);
  // Slopes of the fitted line
  const double tanAlpha = lineDir.x() / lineDir.z();
  const double tanBeta = lineDir.y() / lineDir.z();
  const double rho2 = Acts::square(tanAlpha) + Acts::square(tanBeta);
  // A direction along z has no phi
  if (rho2 < std::numeric_limits<double>::epsilon()) {
    return;
  }
  const double rho = std::sqrt(rho2);
  // Jacobian of (theta, phi) w.r.t. (tanBeta, tanAlpha)
  const double dThdTb = tanBeta / (rho * (1. + rho2));
  const double dThdTa = tanAlpha / (rho * (1. + rho2));
  const double dPhdTb = tanAlpha / rho2;
  const double dPhdTa = -tanBeta / rho2;
  const double varTanB = Acts::square(precResult.dTheta);
  const double varTanA = Acts::square(nonPrecResult.dTheta);
  cov(iY0, iY0) = Acts::square(precResult.dY0);
  cov(iX0, iX0) = Acts::square(nonPrecResult.dY0);
  cov(iTh, iTh) =
      Acts::square(dThdTb) * varTanB + Acts::square(dThdTa) * varTanA;
  cov(iPh, iPh) =
      Acts::square(dPhdTb) * varTanB + Acts::square(dPhdTa) * varTanA;
  cov(iTh, iPh) = cov(iPh, iTh) =
      dThdTb * dPhdTb * varTanB + dThdTa * dPhdTa * varTanA;
  cov(iY0, iTh) = cov(iTh, iY0) = precResult.covY0Theta * dThdTb;
  cov(iY0, iPh) = cov(iPh, iY0) = precResult.covY0Theta * dPhdTb;
  cov(iX0, iTh) = cov(iTh, iX0) = nonPrecResult.covY0Theta * dThdTa;
  cov(iX0, iPh) = cov(iPh, iX0) = nonPrecResult.covY0Theta * dPhdTa;
}

std::vector<CompositeSpacePointLineFitter::FitParIndex>
CompositeSpacePointLineFitter::extractFitablePars(
    const DoFcounts& hitCounts) const {
  if (!m_cfg.parsToUse.empty()) {
    return m_cfg.parsToUse;
  }
  std::vector<FitParIndex> pars{};
  using enum FitParIndex;
  if (hitCounts.nonBending > 1) {
    pars.insert(pars.end(), {x0, phi});
  }
  // Measurements in the bending direction
  if (hitCounts.bending > 1) {
    pars.insert(pars.end(), {y0, theta});
  }
  // Time measurements
  if (m_cfg.fitT0 && (hitCounts.time + hitCounts.straw) > 1) {
    pars.push_back(t0);
  }
  std::ranges::sort(pars);
  return pars;
}

bool CompositeSpacePointLineFitter::withinRange(
    const double parValue, const FitParIndex fitPar) const {
  const auto p = toUnderlying(fitPar);
  /// The parameter is unbound
  if (m_cfg.ranges[p][0] >= m_cfg.ranges[p][1]) {
    return true;
  }
  if (parValue < m_cfg.ranges[p][0] || parValue > m_cfg.ranges[p][1]) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": The parameter "
                          << detail::CompSpacePointAuxiliaries::parName(fitPar)
                          << " " << parValue << " is out range ["
                          << m_cfg.ranges[p][0] << ";" << m_cfg.ranges[p][1]
                          << "]");
    return false;
  }
  return true;
}
bool CompositeSpacePointLineFitter::withinRange(
    const ParamVec_t& pars, const std::vector<FitParIndex>& parsToUse) const {
  return std::ranges::all_of(parsToUse, [&](const FitParIndex par) {
    return withinRange(pars[toUnderlying(par)], par);
  });
}

void CompositeSpacePointLineFitter::FitParameters::print(
    std::ostream& ostr) const {
  using namespace Acts::UnitLiterals;
  using enum FitParIndex;
  ostr << "parameters converged: " << (converged ? "yes" : "no") << ", ";
  ostr << "# iterations: " << nIter << ", ";
  ostr << "chi2: " << chi2 << ", ";
  ostr << "parameters - ";
  ostr << std::format("theta: {:.3f}, ",
                      parameters[toUnderlying(theta)] / 1._degree);
  ostr << std::format("phi: {:.3f}, ",
                      parameters[toUnderlying(phi)] / 1._degree);
  ostr << std::format("y0: {:.3f}, ", parameters[toUnderlying(y0)]);
  ostr << std::format("x0: {:.3f}, ", parameters[toUnderlying(x0)]);
  ostr << std::format("t0: {:.3f}, ", parameters[toUnderlying(t0)] / 1._ns);
  ostr << "covariance - \n" << covariance << ",\n";
}
}  // namespace Acts::Experimental
