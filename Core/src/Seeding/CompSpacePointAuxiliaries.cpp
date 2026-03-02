// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
//
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

namespace Acts::Experimental::detail {

using Vector = CompSpacePointAuxiliaries::Vector;

namespace {
double angle(const Vector& v1, const Vector& v2) {
  using namespace Acts::UnitLiterals;
  return std::acos(std::clamp(v1.dot(v2), -1., 1.)) / 1_degree;
}
constexpr double s_tolerance = 1.e-12;
constexpr auto bendingComp =
    toUnderlying(CompSpacePointAuxiliaries::ResidualIdx::bending);
constexpr auto nonBendingComp =
    toUnderlying(CompSpacePointAuxiliaries::ResidualIdx::nonBending);
constexpr auto timeComp =
    toUnderlying(CompSpacePointAuxiliaries::ResidualIdx::time);
}  // namespace

std::string CompSpacePointAuxiliaries::parName(const FitParIndex idx) {
  switch (idx) {
    using enum FitParIndex;
    case x0:
      return "x0";
    case y0:
      return "y0";
    case theta:
      return "theta";
    case phi:
      return "phi";
    case t0:
      return "t0";
    case nPars:
      return "nPars";
  }
  return "unknown";
}

void CompSpacePointAuxiliaries::ChiSqWithDerivatives::reset() {
  *this = ChiSqWithDerivatives();
}

CompSpacePointAuxiliaries::CompSpacePointAuxiliaries(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {
  std::ranges::sort(m_cfg.parsToUse);
}
void CompSpacePointAuxiliaries::updateChiSq(
    ChiSqWithDerivatives& chiSqObj, const std::array<double, 3>& cov) const {
  // Calculate the inverted covariacne
  std::array<double, 3> invCov{cov};
  for (auto& val : invCov) {
    val = val > s_tolerance ? 1. / val : 0.;
  }

  auto contract = [&invCov, this](const Vector& v1, const Vector& v2) {
    const double term = v1[0] * v2[0] * invCov[0] + v1[1] * v2[1] * invCov[1] +
                        v1[2] * v2[2] * invCov[2];
    ACTS_VERBOSE("updateChiSq() - Contribution from "
                 << toString(v1) << " & " << toString(v2) << " is ("
                 << (v1[0] * v2[0] * invCov[0]) << ", "
                 << (v1[1] * v2[1] * invCov[1]) << ", "
                 << v1[2] * v2[2] * invCov[2] << ") -> overall: " << term);
    return term;
  };
  ACTS_VERBOSE("updateChiSq() - Update chi2 value.");
  chiSqObj.chi2 += contract(residual(), residual());
  for (const auto partial1 : m_cfg.parsToUse) {
    ACTS_VERBOSE("updateChiSq() - Update chi2 derivative w.r.t. "
                 << parName(partial1) << ".");
    chiSqObj.gradient[toUnderlying(partial1)] +=
        2. * contract(residual(), gradient(partial1));
    for (const auto partial2 : m_cfg.parsToUse) {
      if (partial2 > partial1) {
        break;
      }
      ACTS_VERBOSE("updateChiSq() - Update chi2 hessian w.r.t. "
                   << parName(partial1) << " & " << parName(partial2) << ".");
      chiSqObj.hessian(toUnderlying(partial1), toUnderlying(partial2)) +=
          2. * contract(gradient(partial1), gradient(partial2)) +
          (m_cfg.useHessian
               ? 2. * contract(residual(), hessian(partial1, partial2))
               : 0.0);
    }
  }
}
void CompSpacePointAuxiliaries::symmetrizeHessian(
    ChiSqWithDerivatives& chiSqObj) const {
  for (const auto partial1 : m_cfg.parsToUse) {
    for (const auto partial2 : m_cfg.parsToUse) {
      if (partial2 >= partial1) {
        break;
      }
      chiSqObj.hessian(toUnderlying(partial2), toUnderlying(partial1)) =
          chiSqObj.hessian(toUnderlying(partial1), toUnderlying(partial2));
    }
  }
}
const Vector& CompSpacePointAuxiliaries::residual() const {
  return m_residual;
}
const Vector& CompSpacePointAuxiliaries::gradient(
    const FitParIndex param) const {
  assert(toUnderlying(param) < m_gradient.size());
  return m_gradient[toUnderlying(param)];
}
const Vector& CompSpacePointAuxiliaries::hessian(
    const FitParIndex param1, const FitParIndex param2) const {
  const auto idx =
      vecIdxFromSymMat<s_nPars>(toUnderlying(param1), toUnderlying(param2));
  assert(idx < m_hessian.size());
  return m_hessian[idx];
}

void CompSpacePointAuxiliaries::reset() {
  m_residual.setZero();
  for (Vector3& grad : m_gradient) {
    grad.setZero();
  }
  if (m_cfg.useHessian) {
    for (Vector3& hess : m_hessian) {
      hess.setZero();
    }
  }
}
void CompSpacePointAuxiliaries::resetTime() {
  m_gradient[toUnderlying(FitParIndex::t0)].setZero();
  if (m_cfg.useHessian) {
    for (const auto partial : m_cfg.parsToUse) {
      const auto pIdx = vecIdxFromSymMat<s_nPars>(toUnderlying(FitParIndex::t0),
                                                  toUnderlying(partial));
      m_hessian[pIdx].setZero();
    }
  }
}

bool CompSpacePointAuxiliaries::updateStrawResidual(const Line_t& line,
                                                    const Vector& hitMinSeg,
                                                    const Vector& wireDir,
                                                    const double driftRadius) {
  // Fetch the hit position & direction
  if (!updateStrawAuxiliaries(line, wireDir)) {
    ACTS_WARNING("updateStrawResidual() - The line "
                 << toString(line.direction())
                 << " is parallel to the measurement: " << toString(wireDir));
    return false;
  }
  // Distance from the segment line to the tube wire
  const double lineDist = m_projDir.cross(wireDir).dot(hitMinSeg);
  const double resVal = lineDist - driftRadius;
  m_residual = resVal * Vector::Unit(bendingComp);

  // Calculate the first derivative of the residual
  for (const auto partial : m_cfg.parsToUse) {
    const auto pIdx = toUnderlying(partial);
    if (isDirectionParam(partial)) {
      const double partialDist =
          m_gradProjDir[pIdx].cross(wireDir).dot(hitMinSeg);
      ACTS_VERBOSE("updateStrawResidual() - Partial "
                   << parName(partial) << ": " << toString(m_gradProjDir[pIdx])
                   << ", residual grad: " << partialDist);
      m_gradient[pIdx] = partialDist * Vector::Unit(bendingComp);

    } else if (isPositionParam(partial)) {
      const double partialDist = -m_projDir.cross(wireDir).dot(
          line.gradient(static_cast<LineIndex>(partial)));
      ACTS_VERBOSE("updateStrawResidual() - Partial " << parName(partial)
                                                      << ": " << partialDist);
      m_gradient[pIdx] = partialDist * Vector3::Unit(bendingComp);
    }
  }

  if (!m_cfg.useHessian) {
    ACTS_VERBOSE("updateStrawResidual() - Skip Hessian calculation");
    return true;
  }
  // Loop to include the second order derivatvies
  for (const auto partial1 : m_cfg.parsToUse) {
    if (partial1 == FitParIndex::t0) {
      continue;
    }
    for (const auto partial2 : m_cfg.parsToUse) {
      if (partial2 == FitParIndex::t0) {
        continue;
      } else if (partial2 > partial1) {
        break;
      }
      ACTS_VERBOSE("updateStrawResidual() - Calculate Hessian for parameters "
                   << parName(partial1) << ", " << parName(partial2) << ".");
      // Second derivative w.r.t to the position parameters is zero
      const auto hIdx1 = toUnderlying(partial1);
      const auto hIdx2 = toUnderlying(partial2);
      const std::size_t resIdx = vecIdxFromSymMat<s_nPars>(hIdx1, hIdx2);
      if (!(isDirectionParam(partial1) || isDirectionParam(partial2)) ||
          partial2 == FitParIndex::t0 || partial1 == FitParIndex::t0) {
        m_hessian[resIdx].setZero();
        continue;
      }
      const std::size_t projDirHess =
          vecIdxFromSymMat<s_nLinePars>(hIdx1, hIdx2);
      // Pure angular derivatives of the residual
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        const double partialSqDist = m_hessianProjDir[projDirHess].cross(wireDir).dot(hitMinSeg);
        // clang-format on
        m_hessian[resIdx] = partialSqDist * Vector3::Unit(bendingComp);
      } else {
        // Angular & Spatial mixed terms
        const auto angParam = isDirectionParam(partial1) ? hIdx1 : hIdx2;
        const auto posParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial2 : partial1);
        // clang-format off
        m_hessian[resIdx] = -m_gradProjDir[angParam].cross(wireDir).dot(line.gradient(posParam)) * Vector3::Unit(bendingComp);
        // clang-format on
      }
      ACTS_VERBOSE("updateStrawResidual() - Hessian of the residual is "
                   << m_hessian[resIdx][bendingComp]);
    }
  }
  return true;
}
inline bool CompSpacePointAuxiliaries::updateStrawAuxiliaries(
    const Line_t& line, const Vector& wireDir) {
  const Vector& lineDir = line.direction();
  // Between two calls the wire projection has not changed
  const double wireProject = lineDir.dot(wireDir);

  m_wireProject = wireProject;
  const double projDirLenSq = 1. - square(m_wireProject);
  // The line is parallel to the wire
  if (projDirLenSq < s_tolerance) {
    ACTS_VERBOSE("updateStrawAuxiliaries() - Line & wire are parallel: "
                 << toString(wireDir) << " vs. " << toString(lineDir));
    m_invProjDirLenSq = 0.;
    reset();
    return false;
  }
  m_invProjDirLenSq = 1. / projDirLenSq;
  m_invProjDirLen = std::sqrt(m_invProjDirLenSq);
  // Project the segment line onto the wire plane and normalize
  Vector newProjDir = (lineDir - m_wireProject * wireDir) * m_invProjDirLen;
  ACTS_VERBOSE("updateStrawAuxiliaries() - Projected direction is "
               << toString(newProjDir) << ", |D| " << newProjDir.norm());

  if (Acts::abs(newProjDir.dot(m_projDir) - 1.) <
      2. * std::numeric_limits<double>::epsilon()) {
    ACTS_VERBOSE(
        "updateStrawAuxiliaries() - The old & the new dir projection coincide. "
        << "Don't update the auxiliaries");
    return true;
  }

  ACTS_VERBOSE("updateStrawAuxiliaries() - Update variables.");
  m_projDir = std::move(newProjDir);

  // Loop over all configured parameters and calculate the partials
  // of the wire projection
  for (auto partial : m_cfg.parsToUse) {
    // Skip the parameters that are not directional
    if (!isDirectionParam(partial)) {
      continue;
    }
    const std::size_t pIdx = toUnderlying(partial);
    const Vector& lGrad{line.gradient(static_cast<LineIndex>(partial))};
    m_projDirLenPartial[pIdx] = lGrad.dot(wireDir);
    // clang-format off
    m_gradProjDir[pIdx] = m_invProjDirLen * lGrad - (m_invProjDirLen *m_projDirLenPartial[pIdx]) * wireDir +
                          (m_projDirLenPartial[pIdx] * m_wireProject * m_invProjDirLenSq) * m_projDir;
    // clang-format on
    ACTS_VERBOSE("updateStrawAuxiliaries() - First derivative w.r.t. "
                 << parName(partial) << " is "
                 << toString(m_gradProjDir[pIdx]));
    // skip the Hessian calculation, if toggled
    if (!m_cfg.useHessian) {
      continue;
    }
    for (auto partial1 : m_cfg.parsToUse) {
      // Skip the parameters that are not directional or avoid double counting
      if (!isDirectionParam(partial1)) {
        continue;
      } else if (partial1 > partial) {
        break;
      }
      const auto param = toUnderlying(partial);
      const auto param1 = toUnderlying(partial1);
      const auto resIdx = vecIdxFromSymMat<s_nLinePars>(param, param1);
      const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial),
                                            static_cast<LineIndex>(partial1));
      const double partSqLineProject = lHessian.dot(wireDir);

      auto calcMixedTerms = [this](const std::size_t p1, const std::size_t p2) {
        return m_projDirLenPartial[p1] * m_wireProject * m_gradProjDir[p2] +
               0.5 * m_projDirLenPartial[p1] * m_projDirLenPartial[p2] *
                   m_invProjDirLenSq * m_projDir;
      };
      auto calcMixedHessian = [&]() -> Vector {
        if (param1 != param) {
          return calcMixedTerms(param1, param) + calcMixedTerms(param, param1);
        }
        return 2. * calcMixedTerms(param, param1);
      };
      // clang-format off
      m_hessianProjDir[resIdx] = (lHessian - partSqLineProject * wireDir) * m_invProjDirLen +
                               m_invProjDirLenSq * (calcMixedHessian() + (partSqLineProject * m_wireProject) * m_projDir);
      // clang-format on
      ACTS_VERBOSE("updateStrawAuxiliaries() - Hessian w.r.t. "
                   << parName(partial) << ", " << parName(partial1) << " is "
                   << toString(m_hessianProjDir[resIdx]));
    }
  }
  return true;
}

void CompSpacePointAuxiliaries::updateAlongTheStraw(const Line_t& line,
                                                    const Vector& hitMinSeg,
                                                    const Vector& wireDir) {
  // clang-format off
  m_residual[nonBendingComp] = m_invProjDirLenSq *
                            hitMinSeg.dot(m_wireProject * line.direction() - wireDir);
  // clang-format on
  ACTS_VERBOSE("updateAlongTheStraw() - Residual is "
               << m_residual[nonBendingComp]);

  for (const auto partial : m_cfg.parsToUse) {
    if (partial == FitParIndex::t0) {
      continue;
    }
    const auto param = toUnderlying(partial);
    const Vector& lGradient = line.gradient(static_cast<LineIndex>(partial));
    if (isDirectionParam(partial)) {
      // clang-format off
      m_gradient[param][nonBendingComp] = m_invProjDirLenSq * (
                                       hitMinSeg.dot(m_wireProject * lGradient + m_projDirLenPartial[param] * line.direction()) +
                                       2. * m_residual[nonBendingComp] * m_wireProject * m_projDirLenPartial[param]);
      // clang-format on

    } else if (isPositionParam(partial)) {
      // clang-format off
      m_gradient[param][nonBendingComp] = -m_invProjDirLenSq * (lGradient.dot(m_wireProject * line.direction() - wireDir));
      // clang-format on
    }
    ACTS_VERBOSE("updateAlongTheStraw() - First derivative w.r.t "
                 << parName(partial) << " is "
                 << m_gradient[param][nonBendingComp]);
  }
  if (!m_cfg.useHessian) {
    return;
  }

  for (auto partial1 : m_cfg.parsToUse) {
    if (partial1 == FitParIndex::t0) {
      continue;
    }
    for (auto partial2 : m_cfg.parsToUse) {
      if (partial2 == FitParIndex::t0) {
        continue;
      } else if (partial2 > partial1) {
        break;
      }
      // second derivative w.r.t intercepts is always 0
      if (!(isDirectionParam(partial2) || isDirectionParam(partial1)) ||
          partial2 == FitParIndex::t0 || partial1 == FitParIndex::t0) {
        continue;
      }
      const auto param = toUnderlying(partial1);
      const auto param1 = toUnderlying(partial2);
      const auto hessIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        auto calcMixedHessian = [this, &hitMinSeg, &line](const std::size_t p1,
                                                          const std::size_t p2) -> double {
          const double term1 = hitMinSeg.dot(line.gradient(static_cast<LineIndex>(p1))) * m_projDirLenPartial[p2];
          const double term2 = m_residual[nonBendingComp] * m_projDirLenPartial[p1] * m_projDirLenPartial[p2];
          const double term3 = 2. * m_wireProject * m_projDirLenPartial[p1] * m_gradient[p2][nonBendingComp];
          return (term1 + term2 + term3) * m_invProjDirLenSq;
        };
        const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial1), static_cast<LineIndex>(partial2));
        const double projHess = lHessian.dot(wireDir);
        m_hessian[hessIdx][nonBendingComp] = (param != param1 ? calcMixedHessian(param, param1) + calcMixedHessian(param1, param)
                                                              : 2. * calcMixedHessian(param1, param)) +
                                        m_invProjDirLenSq * (hitMinSeg.dot(m_wireProject * lHessian  + projHess* line.direction()) +
                                                             2. * m_residual[nonBendingComp] * m_wireProject * projHess);
        // clang-format on
      } else {
        // Mixed positional and angular derivative
        const auto angParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial1 : partial2);
        const auto posParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial2 : partial1);
        const auto angIdx = toUnderlying(angParam);
        // clang-format off
        m_hessian[hessIdx][nonBendingComp] = m_invProjDirLenSq * (
                  2.* (m_gradient[toUnderlying(posParam)][nonBendingComp] * m_wireProject * m_projDirLenPartial[angIdx])
                  - line.gradient(posParam).dot(m_wireProject * line.gradient(angParam) + m_projDirLenPartial[angIdx] * line.direction()));
        // clang-format on
      }
      ACTS_VERBOSE("updateAlongTheStraw() - Second derivative w.r.t. "
                   << parName(partial1) << ", " << parName(partial2) << " is. "
                   << m_hessian[hessIdx][nonBendingComp]);
    }
  }
}

void CompSpacePointAuxiliaries::updateStripResidual(
    const Line_t& line, const Vector& normal, const Vector& sensorN,
    const Vector& sensorD, const Vector& stripPos, const bool isBending,
    const bool isNonBending) {
  const double normDot = normal.dot(line.direction());

  constexpr double tolerance = 1.e-12;
  if (std::abs(normDot) < tolerance) {
    reset();
    ACTS_WARNING(
        "updateStripResidual() - Segment line is embedded into the strip plane "
        << toString(line.direction()) << ", normal: " << toString(normal));
    return;
  }
  const double planeOffSet = normal.dot(stripPos);
  ACTS_VERBOSE("updateStripResidual() - Plane normal: "
               << toString(normal) << " |" << normal.norm()
               << "|, line direction: " << toString(line.direction())
               << ", angle: " << angle(normal, line.direction()));
  const double invNormDot = 1. / normDot;
  const double travelledDist =
      (planeOffSet - line.position().dot(normal)) * invNormDot;
  ACTS_VERBOSE("updateStripResidual() - Need to travel "
               << travelledDist << " mm to reach the strip plane. ");

  const Vector& b1 = isBending ? sensorN : sensorD;
  const Vector& b2 = isBending ? sensorD : sensorN;

  const double b1DotB2 = b1.dot(b2);
  if (Acts::abs(Acts::abs(b1DotB2) - 1.) < tolerance) {
    ACTS_WARNING("updateStripResidual() The two vectors "
                 << toString(b1) << ", " << toString(b2) << " with angle: "
                 << angle(b1, b2) << " don't span the plane.");
    reset();
    return;
  }
  // Linear independent vectors
  // Normal vector is indeed normal onto the plane spaned by these two vectors
  assert(Acts::abs(b1.dot(normal)) < tolerance);
  assert(Acts::abs(b2.dot(normal)) < tolerance);
  assert(isBending || isNonBending);
  const bool decompStereo =
      (m_cfg.calcAlongStrip || (isBending && isNonBending)) &&
      Acts::abs(b1DotB2) > tolerance;
  if (decompStereo) {
    // A = lambda B1 + kappa B2
    //                <A, B2> - <A,B1><B1,B2>
    //   --> kappa = --------------------------
    //                    1 - <B1,B2>^{2}
    //                <A, B1> - <A,B2><B1,B2>
    //   --> lambda = --------------------------
    //                    1 - <B1,B2>^{2}
    const double invDist = 1. / (1. - square(b1DotB2));
    m_stereoTrf(bendingComp, bendingComp) =
        m_stereoTrf(nonBendingComp, nonBendingComp) = invDist;
    m_stereoTrf(bendingComp, nonBendingComp) =
        m_stereoTrf(nonBendingComp, bendingComp) = -b1DotB2 * invDist;
  }
  ACTS_VERBOSE("updateStripResidual() - Measures bending "
               << (isBending ? "yes" : "no") << ", measures non-bending "
               << (isNonBending ? "yes" : "no"));
  ACTS_VERBOSE("updateStripResidual() - Strip orientation b1: "
               << toString(b1) << " |" << b1.norm() << "|"
               << ", b2: " << toString(b2) << " |" << b2.norm() << "|"
               << ", stereo angle: " << angle(b1, b2));

  auto assignResidual = [&](const Vector& calcDistance, Vector& residual) {
    residual[bendingComp] =
        isBending || m_cfg.calcAlongStrip ? b1.dot(calcDistance) : 0.;
    residual[nonBendingComp] =
        isNonBending || m_cfg.calcAlongStrip ? b2.dot(calcDistance) : 0.;
    if (decompStereo) {
      auto spatial = residual.block<2, 1>(0, 0);
      spatial = m_stereoTrf * spatial;
    }
    residual[timeComp] = 0.;
    ACTS_VERBOSE(
        "updateStripResidual() - Distance: "
        << toString(calcDistance) << ", <calc, n> =" << calcDistance.dot(normal)
        << " -> residual: " << toString(residual) << " --> closure test:"
        << toString(residual[bendingComp] * b1 +
                    residual[nonBendingComp] * b2));
  };
  // Update the residual accordingly
  assignResidual(line.position() + travelledDist * line.direction() - stripPos,
                 m_residual);
  for (const auto partial : m_cfg.parsToUse) {
    ACTS_VERBOSE("updateStripResidual() - Calculate partial derivative w.r.t "
                 << parName(partial));
    switch (partial) {
      case FitParIndex::phi:
      case FitParIndex::theta: {
        // clang-format off
        const Vector& lGrad = line.gradient(static_cast<LineIndex>(partial));
        const double partialDist = -travelledDist * invNormDot * normal.dot(lGrad);
        const Vector gradCmp = travelledDist * lGrad +
                               partialDist * line.direction();
        // clang-format on
        assignResidual(gradCmp, m_gradient[toUnderlying(partial)]);
        break;
      }
      case FitParIndex::y0:
      case FitParIndex::x0: {
        // clang-format off
        const Vector& lGrad = line.gradient(static_cast<LineIndex>(partial));
        const Vector gradCmp = lGrad - lGrad.dot(normal) * invNormDot * line.direction();
        assignResidual(gradCmp, m_gradient[toUnderlying(partial)]);
        // clang-format on
        break;
      }
      default:
        // Don't calculate anything for time
        break;
    }
  }

  if (!m_cfg.useHessian) {
    return;
  }
  for (const auto partial1 : m_cfg.parsToUse) {
    for (const auto partial2 : m_cfg.parsToUse) {
      // Avoid double counting
      if (partial2 > partial1) {
        break;
      }
      const auto param = toUnderlying(partial1);
      const auto param1 = toUnderlying(partial2);
      const auto resIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      // At least one parameter needs to be a directional one
      if (!(isDirectionParam(partial1) || isDirectionParam(partial2)) ||
          partial2 == FitParIndex::t0 || partial1 == FitParIndex::t0) {
        m_hessian[resIdx].setZero();
        continue;
      }
      ACTS_VERBOSE(
          "updateStripResidual() - Calculate second partial derivative w.r.t "
          << parName(partial1) << ", " << parName(partial2));

      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        auto calcMixedTerms = [&line, &normal, &invNormDot, &b1, &b2, this](const std::size_t p1,
                                                                            const std::size_t p2) {

          return -normal.dot(line.gradient(static_cast<LineIndex>(p1))) * invNormDot *
                          (m_gradient[p2][bendingComp] * b1 + m_gradient[p2][nonBendingComp] * b2);
        };
        const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial1),
                                              static_cast<LineIndex>(partial2));
        const Vector hessianCmp = travelledDist * (lHessian - normal.dot(lHessian) * invNormDot * line.direction()) +
                                  calcMixedTerms(param1, param) + calcMixedTerms(param, param1);
        assignResidual(hessianCmp, m_hessian[resIdx]);
        // clang-format on
      } else {
        // clang-format off
        const auto angParam = static_cast<LineIndex>(isDirectionParam(partial1) ? partial1 : partial2);
        const auto posParam = static_cast<LineIndex>(isDirectionParam(partial1) ? partial2 : partial1);
        const double lH = (normal.dot(line.gradient(posParam)) * invNormDot);

        const Vector hessianCmp = lH * ((normal.dot(line.gradient(angParam)) * invNormDot) * line.direction() -
                                        line.gradient(angParam));
        // clang-format on
        assignResidual(hessianCmp, m_hessian[resIdx]);
      }
    }
  }
}

void CompSpacePointAuxiliaries::updateTimeStripRes(
    const Vector& sensorN, const Vector& sensorD, const Vector& stripPos,
    const bool isBending, const double recordTime, const double timeOffset) {
  const Vector& b1 = isBending ? sensorN : sensorD;
  const Vector& b2 = isBending ? sensorD : sensorN;

  auto positionInPlane = [&b1, &b2](const Vector& residVec) {
    return b1 * residVec[bendingComp] + b2 * residVec[nonBendingComp];
  };

  // Calculate the line intersection in the global frame
  const Vector globIsect =
      m_cfg.includeToF
          ? m_cfg.localToGlobal * (positionInPlane(residual()) + stripPos)
          : Vector::Zero();
  // To calculate the time of flight
  const double globDist = globIsect.norm();
  const double ToF = globDist / PhysicalConstants::c + timeOffset;
  ACTS_VERBOSE("Global intersection point: "
               << toString(globIsect) << " -> distance: " << globDist
               << " -> time of flight takinig t0 into account: " << ToF);

  m_residual[timeComp] = recordTime - ToF;
  constexpr auto timeIdx = toUnderlying(FitParIndex::t0);

  const double invDist = m_cfg.includeToF && globDist > s_tolerance
                             ? 1. / (PhysicalConstants::c * globDist)
                             : 0.;
  for (const auto partial1 : m_cfg.parsToUse) {
    if (partial1 == FitParIndex::t0) {
      m_gradient[timeIdx] = -Vector::Unit(timeComp);
    }
    // Time component of the spatial residual needs to be updated
    else if (m_cfg.includeToF) {
      Vector& gradVec = m_gradient[toUnderlying(partial1)];
      gradVec[timeComp] = -globIsect.dot(m_cfg.localToGlobal.linear() *
                                         positionInPlane(gradVec)) *
                          invDist;
      ACTS_VERBOSE("Partial of the time residual  w.r.t. "
                   << parName(partial1) << ": " << gradVec[timeComp] << ".");
    }
  }

  if (!m_cfg.useHessian || !m_cfg.includeToF) {
    return;
  }
  for (const auto partial1 : m_cfg.parsToUse) {
    for (const auto partial2 : m_cfg.parsToUse) {
      if (partial2 > partial1) {
        break;
      }
      const auto param1 = toUnderlying(partial1);
      const auto param2 = toUnderlying(partial2);
      Vector& hessVec = m_hessian[vecIdxFromSymMat<s_nPars>(param1, param2)];
      if (partial1 != FitParIndex::t0) {
        // clang-format off
        hessVec[timeComp] = -( globIsect.dot(m_cfg.localToGlobal.linear()*positionInPlane(hessVec))  +
                           positionInPlane(gradient(partial1)).dot(positionInPlane(gradient(partial2)))) * invDist
                        + m_gradient[param1][timeComp] * m_gradient[param2][timeComp] * invDist;
        // clang-format on
        ACTS_VERBOSE("Hessian of the time residual w.r.t. "
                     << parName(partial1) << ", " << parName(partial2) << ": "
                     << hessVec[timeComp] << ".");
      } else {
        hessVec.setZero();
      }
    }
  }
}

void CompSpacePointAuxiliaries::updateTimeStrawRes(
    const Line_t& line, const Vector& hitMinSeg, const Vector& wireDir,
    const double driftR, const double driftV, const double driftA) {
  using namespace Acts::detail::LineHelper;
  using namespace Acts::UnitLiterals;

  const double dSign = std::copysign(1., driftR);
  // Only assign drift velocity and acceleration
  if (!m_cfg.includeToF) {
    resetTime();
    constexpr auto timeIdx = toUnderlying(FitParIndex::t0);
    constexpr auto hessIdx = vecIdxFromSymMat<s_nPars>(timeIdx, timeIdx);
    m_gradient[timeIdx] = dSign * driftV * Vector::Unit(bendingComp);
    m_hessian[hessIdx] = -dSign * driftA * Vector::Unit(bendingComp);
    return;
  }

  const double travDist = hitMinSeg.dot(m_projDir * m_invProjDirLen);

  const Vector clApproach = line.point(travDist);
  const Vector globApproach = m_cfg.localToGlobal * clApproach;
  ACTS_VERBOSE("updateTimeStrawRes() - Point of closest approach along line "
               << toString(clApproach));

  const double distance = globApproach.norm();
  const double ToF = distance / PhysicalConstants::c;
  ACTS_VERBOSE("updateTimeStrawRes() - Distance from the global origin: "
               << distance << " -> time of flight: " << ToF / 1._ns);
  ACTS_VERBOSE("updateTimeStrawRes() - Drift radius: "
               << driftR << ", driftV: " << driftV << ", driftA: " << driftA);

  const double invDist =
      distance > s_tolerance ? 1. / (distance * PhysicalConstants::c) : 0.;
  for (const auto partial : m_cfg.parsToUse) {
    const auto idx = toUnderlying(partial);
    if (partial == FitParIndex::t0) {
      m_gradient[idx] = dSign * driftV * Vector::Unit(bendingComp);
      continue;
    }
    const Vector3& lGrad = line.gradient(static_cast<LineIndex>(partial));
    if (isPositionParam(partial)) {
      // clang-format off
      m_gradCloseApproach[idx] = lGrad - lGrad.dot(m_projDir * m_invProjDirLen) * line.direction();
      // clang-format on
    } else if (isDirectionParam(partial)) {
      // clang-format off
      m_gradCloseApproach[idx] =  hitMinSeg.dot(m_gradProjDir[idx]) * m_invProjDirLen * line.direction()
                        + travDist * lGrad
                        + m_wireProject  * m_projDirLenPartial[idx] * m_invProjDirLenSq * travDist  * line.direction();
      // clang-format on
    }
    m_partialApproachDist[idx] = globApproach.dot(m_cfg.localToGlobal.linear() *
                                                  m_gradCloseApproach[idx]) *
                                 invDist;
    ACTS_VERBOSE("updateTimeStrawRes() - Correct the partial derivative w.r.t. "
                 << parName(partial) << " " << m_gradient[idx][bendingComp]
                 << " by " << toString(m_gradCloseApproach[idx])
                 << " dCoA: " << m_partialApproachDist[idx] * driftV << " -> "
                 << (m_gradient[idx][bendingComp] +
                     dSign * m_partialApproachDist[idx] * driftV));
    m_gradient[idx][bendingComp] -=
        -dSign * m_partialApproachDist[idx] * driftV;
  }

  if (!m_cfg.useHessian) {
    return;
  }

  for (const auto partial1 : m_cfg.parsToUse) {
    for (const auto partial2 : m_cfg.parsToUse) {
      if (partial2 > partial1) {
        break;
      }
      const auto idx1 = toUnderlying(partial1);
      const auto idx2 = toUnderlying(partial2);
      const auto hessIdx = vecIdxFromSymMat<s_nPars>(idx1, idx2);
      if (partial1 == FitParIndex::t0) {
        if (partial2 == FitParIndex::t0) {
          m_hessian[hessIdx] = -dSign * driftA * Vector::Unit(bendingComp);
        } else {
          m_hessian[hessIdx] = -dSign * driftA * m_partialApproachDist[idx2] *
                               Vector::Unit(bendingComp);
        }
        ACTS_VERBOSE("updateTimeStrawRes() -"
                     << " Second partial derivative of the drift "
                        "radius w.r.t. "
                     << parName(partial1) << ", " << parName(partial2) << " is "
                     << toString(m_hessian[hessIdx]));
        continue;
      }
      Vector hessCmp{Vector::Zero()};
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        auto calcMixedTerms = [this, &hitMinSeg, &travDist, &line]
           (const std::size_t pIdx1, const std::size_t pIdx2) ->Vector {
           return m_projDirLenPartial[pIdx1] * m_wireProject * m_invProjDirLenSq * m_gradCloseApproach[pIdx2]
                  + hitMinSeg.dot(m_gradProjDir[pIdx1]) * m_invProjDirLen * line.gradient(static_cast<LineIndex>(pIdx2))
                  + 0.5 * m_invProjDirLenSq *  m_projDirLenPartial[pIdx1] * m_projDirLenPartial[pIdx2] *
                   travDist*( 1. + Acts::pow(m_wireProject* m_invProjDirLen, 2)) * line.direction();
        };

        const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial1),
                                              static_cast<LineIndex>(partial2));

        const std::size_t hessProjIdx = vecIdxFromSymMat<s_nLinePars>(idx1, idx2);
        hessCmp = m_invProjDirLen * hitMinSeg.dot(m_hessianProjDir[hessProjIdx]) * line.direction()
                 + travDist*lHessian
                 + travDist*m_wireProject*lHessian.dot(wireDir)*m_invProjDirLenSq *line.direction();
        // clang-format on
        if (idx1 != idx2) {
          hessCmp += calcMixedTerms(idx1, idx2) + calcMixedTerms(idx2, idx1);
        } else {
          hessCmp += 2. * calcMixedTerms(idx1, idx2);
        }
      } else if (isDirectionParam(partial1) || isDirectionParam(partial2)) {
        // clang-format off
        const auto angParam = isDirectionParam(partial1) ? idx1 : idx2;
        const Vector& posPartial = line.gradient(static_cast<LineIndex>(isDirectionParam(partial1) ? partial2 : partial1));
        const Vector& angPartial = line.gradient(static_cast<LineIndex>(angParam));
        const double gradProjDist = - posPartial.dot(m_projDir) * m_invProjDirLen;
        hessCmp = - posPartial.dot(m_gradProjDir[angParam]) * m_invProjDirLen * line.direction()
                + gradProjDist * angPartial
                + m_wireProject  * m_projDirLenPartial[angParam] *m_invProjDirLenSq * gradProjDist * line.direction();
        // clang-format on
      }
      // clang-format off
      const double partialCoA = m_gradCloseApproach[idx1].dot(m_gradCloseApproach[idx2]) * invDist +
                                globApproach.dot(m_cfg.localToGlobal.linear()*hessCmp) * invDist -
                                m_partialApproachDist[idx1] * m_partialApproachDist[idx2]* invDist * Acts::pow(PhysicalConstants::c,2);
      const double hessianR = driftA * m_partialApproachDist[idx1] * m_partialApproachDist[idx2]
                            - driftV * partialCoA;
      // clang-format on
      ACTS_VERBOSE(
          "updateTimeStrawRes() - Correct the second partial derivative w.r.t "
          << parName(partial1) << ", " << parName(partial2) << " "
          << m_hessian[hessIdx][bendingComp] << " by hessianR: " << hessianR
          << ", partialCoA: " << partialCoA << " -> "
          << m_hessian[hessIdx][bendingComp] - dSign * hessianR);

      m_hessian[hessIdx][bendingComp] -= dSign * hessianR;
    }
  }
}

}  // namespace Acts::Experimental::detail
