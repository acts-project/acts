// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

namespace Acts::Experimental {

CompositeSpacePointLineFitter::CompositeSpacePointLineFitter(
    const Config& cfg, std::unique_ptr<const Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {}

detail::FastStrawLineFitter::Config
CompositeSpacePointLineFitter::fastFitterCfg() const {
  detail::FastStrawLineFitter::Config fastCfg{};
  fastCfg.maxIter = m_cfg.maxIter;
  fastCfg.precCutOff = m_cfg.precCutOff;
  return fastCfg;
}

std::vector<CompositeSpacePointLineFitter::FitParIndex>
CompositeSpacePointLineFitter::extractFitablePars(
    const std::array<std::size_t, 3>& hitCounts) {
  std::vector<FitParIndex> pars{};
  using enum FitParIndex;
  const auto& [nLoc0, nLoc1, nTime] = hitCounts;
  if (nLoc0 > 1) {
    pars.insert(pars.end(), {x0, phi});
  }
  // Measurements in the bending direction
  if (nLoc1 > 1) {
    pars.insert(pars.end(), {y0, theta});
  }
  // Time measurements
  if (nTime > 1) {
    pars.push_back(t0);
  }
  std::ranges::sort(pars);
  return pars;
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
