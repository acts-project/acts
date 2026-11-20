// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include <optional>

#include <Eigen/Cholesky>

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

bool CompositeSpacePointLineFitter::evaluateCovariance(
    const CovMat_t& gaussNewtonHessian, std::vector<FitParIndex> parsToUse,
    CovMat_t& covariance) const {
  using DynMat_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0,
                                 s_nPars, s_nPars>;
  /// Returns 2 H^-1 of the Hessian block of the given parameters, or
  /// std::nullopt if the block is not positive definite. The block is
  /// normalized to unit diagonal before the inversion, as its entries span
  /// many orders of magnitude (the time entries are tiny compared to the
  /// angular ones)
  auto invertBlock =
      [&gaussNewtonHessian](
          const std::vector<FitParIndex>& pars) -> std::optional<DynMat_t> {
    DynMat_t block(pars.size(), pars.size());
    for (std::size_t i = 0; i < pars.size(); ++i) {
      for (std::size_t j = 0; j < pars.size(); ++j) {
        block(i, j) =
            gaussNewtonHessian(toUnderlying(pars[i]), toUnderlying(pars[j]));
      }
    }
    if ((block.diagonal().array() <= 0.).any()) {
      return std::nullopt;
    }
    const auto scale = block.diagonal().cwiseSqrt().cwiseInverse().eval();
    block = scale.asDiagonal() * block * scale.asDiagonal();
    Eigen::LLT<DynMat_t> llt{block};
    if (llt.info() != Eigen::Success) {
      return std::nullopt;
    }
    return DynMat_t{2. * scale.asDiagonal() *
                    llt.solve(DynMat_t::Identity(pars.size(), pars.size())) *
                    scale.asDiagonal()};
  };
  constexpr auto t0Idx = toUnderlying(FitParIndex::t0);
  std::optional<DynMat_t> inverse = invertBlock(parsToUse);
  const bool t0Unconstrained =
      !parsToUse.empty() && parsToUse.back() == FitParIndex::t0 &&
      (!inverse || (*inverse)(parsToUse.size() - 1, parsToUse.size() - 1) >
                       m_cfg.unconstrainedT0Variance);
  if (t0Unconstrained) {
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << ": t0 not constrained by the measurements");
    parsToUse.pop_back();
    inverse = invertBlock(parsToUse);
  }
  if (!inverse) {
    return false;
  }
  for (std::size_t i = 0; i < parsToUse.size(); ++i) {
    for (std::size_t j = 0; j < parsToUse.size(); ++j) {
      covariance(toUnderlying(parsToUse[i]), toUnderlying(parsToUse[j])) =
          (*inverse)(i, j);
    }
  }
  if (t0Unconstrained) {
    for (const auto par : parsToUse) {
      covariance(toUnderlying(par), t0Idx) = 0.;
      covariance(t0Idx, toUnderlying(par)) = 0.;
    }
    covariance(t0Idx, t0Idx) = m_cfg.unconstrainedT0Variance;
  }
  return true;
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
