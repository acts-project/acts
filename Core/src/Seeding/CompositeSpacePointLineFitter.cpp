// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

namespace Acts::Experimental {

namespace {
/// @brief Covariance of a fast fit in its (intercept, angle) parameters
/// @param result: Fast fit result
SquareMatrix<2> planeCov(const detail::FastStrawLineFitter::FitResult& result) {
  SquareMatrix<2> cov = SquareMatrix<2>::Zero();
  cov(0, 0) = Acts::square(result.dY0);
  cov(1, 1) = Acts::square(result.dTheta);
  cov(0, 1) = cov(1, 0) = result.covY0Theta;
  return cov;
}
}  // namespace

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

void CompositeSpacePointLineFitter::fastFitCovariance(
    const FastFitter_t::FitResult& planeResult,
    const FastFitter_t::ResidualIdx projection, const Vector_t& lineDir,
    CovMat_t& cov) {
  using enum FitParIndex;
  assert(projection == FastFitter_t::ResidualIdx::bending ||
         projection == FastFitter_t::ResidualIdx::nonBending);
  const bool bending = projection == FastFitter_t::ResidualIdx::bending;
  const auto iPos = toUnderlying(bending ? y0 : x0);
  const auto iAng = toUnderlying(bending ? theta : phi);
  // The other angle is fixed: dtanBeta / dtheta at fixed phi (bending),
  // dtanAlpha / dphi at fixed theta (non-bending)
  const double dSlopedAng = Line_t::lineToSlopesJacobian(lineDir)(iAng, iAng);
  // The slope does not change with the angle, e.g. tanAlpha for tanBeta = 0
  if (std::abs(dSlopedAng) < std::numeric_limits<double>::epsilon()) {
    return;
  }

  SquareMatrix<2> jac = SquareMatrix<2>::Identity();
  jac(1, 1) = 1. / dSlopedAng;
  cov.block<2, 2>(iPos, iPos) = jac * planeCov(planeResult) * jac.transpose();
}

void CompositeSpacePointLineFitter::fastFitCovariance(
    const FastFitter_t::FitResult& precResult,
    const FastFitter_t::FitResult& nonPrecResult, const Vector_t& lineDir,
    CovMat_t& cov) {
  using enum FitParIndex;
  const SquareMatrix<4> toSlopes = Line_t::lineToSlopesJacobian(lineDir);
  // A direction along z has no phi
  assert(std::abs(toSlopes.determinant()) >
         std::numeric_limits<double>::epsilon());
  const SquareMatrix<4> jac = toSlopes.inverse();

  SquareMatrix<4> planeCovs = SquareMatrix<4>::Zero();
  planeCovs.block<2, 2>(toUnderlying(y0), toUnderlying(y0)) =
      planeCov(precResult);
  planeCovs.block<2, 2>(toUnderlying(x0), toUnderlying(x0)) =
      planeCov(nonPrecResult);
  cov.block<4, 4>(0, 0) = jac * planeCovs * jac.transpose();
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
