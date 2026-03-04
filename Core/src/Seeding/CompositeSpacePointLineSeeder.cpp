// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <format>

namespace {
/// @brief Express an angle in degree
/// @param rad: Angle in radians
constexpr double inDeg(const double rad) {
  using namespace Acts::UnitLiterals;
  return rad / 1._degree;
}
}  // namespace

namespace Acts::Experimental {

using SeedParam_t = CompositeSpacePointLineSeeder::SeedParam_t;
using TangentAmbi = CompositeSpacePointLineSeeder::TangentAmbi;

CompositeSpacePointLineSeeder::CompositeSpacePointLineSeeder(
    const Config& cfg, std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

TangentAmbi CompositeSpacePointLineSeeder::encodeAmbiguity(
    const int signTop, const int signBottom) {
  assert(Acts::abs(signTop) == 1 && Acts::abs(signBottom) == 1);
  using enum TangentAmbi;
  if (signTop == 1 && signBottom == 1) {
    static_assert(s_signCombo[toUnderlying(RR)][0] == 1 &&
                  s_signCombo[toUnderlying(RR)][1] == 1);
    return RR;
  } else if (signTop == 1 && signBottom == -1) {
    static_assert(s_signCombo[toUnderlying(RL)][0] == 1 &&
                  s_signCombo[toUnderlying(RL)][1] == -1);

    return RL;
  } else if (signTop == -1 && signBottom == 1) {
    static_assert(s_signCombo[toUnderlying(LR)][0] == -1 &&
                  s_signCombo[toUnderlying(LR)][1] == 1);
    return LR;
  }
  static_assert(s_signCombo[toUnderlying(LL)][0] == -1 &&
                s_signCombo[toUnderlying(LL)][1] == -1);
  return LL;
}

std::string CompositeSpacePointLineSeeder::toString(const TangentAmbi ambi) {
  switch (ambi) {
    using enum TangentAmbi;
    case RR:
      return "Right - Right";
    case RL:
      return "Right - Left";
    case LR:
      return "Left - Right";
    case LL:
      return "Left - Left";
  }
  return "Undefined";
}

CompositeSpacePointLineSeeder::Line_t CompositeSpacePointLineSeeder::makeLine(
    const SeedParam_t& pars) const {
  using enum ParIdx;
  using Aux_t = detail::CompSpacePointAuxiliaries;
  ACTS_VERBOSE(__func__ << "() - Create line from parameters "
                        << std::format("{:}: {:.2f}, ", Aux_t::parName(x0),
                                       pars[toUnderlying(x0)])
                        << std::format("{:}: {:.2f}, ", Aux_t::parName(y0),
                                       pars[toUnderlying(y0)])
                        << std::format("{:}: {:.3f}, ", Aux_t::parName(theta),
                                       inDeg(pars[toUnderlying(theta)]))
                        << std::format("{:}: {:.3f}", Aux_t::parName(phi),
                                       inDeg(pars[toUnderlying(phi)])));
  return std::make_pair(
      Vector3{pars[toUnderlying(x0)], pars[toUnderlying(y0)], 0.},
      makeDirectionFromPhiTheta(pars[toUnderlying(phi)],
                                pars[toUnderlying(theta)]));
}

SeedParam_t CompositeSpacePointLineSeeder::combineWithPattern(
    const Line_t& tangentSeed, const SeedParam_t& patternParams) const {
  SeedParam_t result{patternParams};
  using enum ParIdx;
  const auto& [seedPos, seedDir] = tangentSeed;
  result[toUnderlying(y0)] = seedPos.y();
  const Vector patternDir = makeLine(result).second;
  const double tanAlpha = patternDir.x() / patternDir.z();
  const double tanBeta = seedDir.y() / seedDir.z();
  const Vector3 dir = makeDirectionFromAxisTangents(tanAlpha, tanBeta);
  result[toUnderlying(phi)] = VectorHelpers::phi(dir);
  result[toUnderlying(theta)] = VectorHelpers::theta(dir);
  return result;
}

void CompositeSpacePointLineSeeder::TwoCircleTangentPars::print(
    std::ostream& ostr) const {
  ostr << toString(ambi);
  ostr << ", "
       << std::format("theta: {:.3f} pm {:.3f}", inDeg(theta), inDeg(dTheta));
  ostr << ", " << std::format("y_{{0}}: {:.2f} : pm {:.2f}", y0, dY0);
}
bool CompositeSpacePointLineSeeder::isValidLine(
    const TwoCircleTangentPars& tangentPars) const {
  if (m_cfg.thetaRange[0] < m_cfg.thetaRange[1] &&
      (tangentPars.theta < m_cfg.thetaRange[0] ||
       tangentPars.theta > m_cfg.thetaRange[1])) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << " - Theta parameter out of range: "
                          << std::format("{:.3f} ,allowed: [{:.3f};{:.3f}].",
                                         inDeg(tangentPars.theta),
                                         inDeg(m_cfg.thetaRange[0]),
                                         inDeg(m_cfg.thetaRange[1])));
    return false;
  }
  if (m_cfg.interceptRange[0] < m_cfg.interceptRange[1] &&
      (tangentPars.y0 < m_cfg.interceptRange[0] ||
       tangentPars.y0 > m_cfg.interceptRange[1])) {
    ACTS_VERBOSE(__func__ << "() " << __LINE__
                          << " - Interscept parameter out of range: "
                          << std::format("{:.2f}, allowed: [{:.2f};{:.2f}].",
                                         tangentPars.y0,
                                         m_cfg.interceptRange[0],
                                         m_cfg.interceptRange[1]));
    return false;
  }
  return true;
}
}  // namespace Acts::Experimental
