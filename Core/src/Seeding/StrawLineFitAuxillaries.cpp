// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
namespace Acts::Experimental::detail {

using Vector = StrawLineFitAuxiliaries::Vector;

namespace {
double angle(const Vector& v1, const Vector& v2) {
  using namespace Acts::UnitLiterals;
  return std::acos(std::clamp(v1.dot(v2), -1., 1.)) / 1_degree;
}
}  // namespace

std::string StrawLineFitAuxiliaries::parName(const FitParIndex idx) {
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
StrawLineFitAuxiliaries::StrawLineFitAuxiliaries(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {
  std::ranges::sort(m_cfg.parsToUse);
}
const Vector& StrawLineFitAuxiliaries::residual() const {
  return m_residual;
}
const Vector& StrawLineFitAuxiliaries::gradient(const FitParIndex param) const {
  const std::uint8_t par = static_cast<std::uint8_t>(param);
  assert(par < m_gradient.size());
  return m_gradient[par];
}
const Vector& StrawLineFitAuxiliaries::hessian(const FitParIndex param1,
                                               const FitParIndex param2) const {
  const std::uint8_t idx = vecIdxFromSymMat<s_nPars>(
      static_cast<std::size_t>(param1), static_cast<std::size_t>(param2));
  assert(idx < m_hessian.size());
  return m_hessian[idx];
}

void StrawLineFitAuxiliaries::reset() {
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

bool StrawLineFitAuxiliaries::updateStrawResidual(const Line_t& line,
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
  /** Distance from the segment line to the tube wire */
  const double lineDist = m_projDir.cross(wireDir).dot(hitMinSeg);
  const double resVal = lineDist - driftRadius;
  m_residual = resVal * Vector::Unit(bending);

  /** Calculate the first derivative of the residual */
  for (const auto partial : m_cfg.parsToUse) {
    const auto pIdx = static_cast<std::uint8_t>(partial);
    if (isDirectionParam(partial)) {
      const double partialDist =
          m_gradProjDir[pIdx].cross(wireDir).dot(hitMinSeg);
      ACTS_VERBOSE("updateStrawResidual() - Partial "
                   << parName(partial) << ": " << toString(m_gradProjDir[pIdx])
                   << ", residual grad: " << partialDist);
      m_gradient[pIdx] = partialDist * Vector::Unit(bending);

    } else if (isPositionParam(partial)) {
      const double partialDist = -m_projDir.cross(wireDir).dot(
          line.gradient(static_cast<LineIndex>(partial)));
      ACTS_VERBOSE("updateStrawResidual() - Partial " << parName(partial)
                                                      << ": " << partialDist);
      m_gradient[pIdx] = partialDist * Vector3::Unit(bending);
    }
  }

  if (!m_cfg.useHessian) {
    ACTS_VERBOSE("updateStrawResidual() - Skip Hessian calculation");
    return true;
  }
  /** Loop to include the second order derivatvies */
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
      if (!(isDirectionParam(partial1) || isDirectionParam(partial2))) {
        continue;
      }
      const auto hIdx1 = static_cast<std::size_t>(partial1);
      const auto hIdx2 = static_cast<std::size_t>(partial2);

      const std::size_t resIdx = vecIdxFromSymMat<s_nLinePars>(hIdx1, hIdx2);
      /// Pure angular derivatives of the residual
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        const double partialSqDist = m_hessianProjDir[resIdx].cross(wireDir).dot(hitMinSeg);
        // clang-format on
        m_hessian[resIdx] = partialSqDist * Vector3::Unit(bending);
      } else {
        /// Angular & Spatial mixed terms
        const auto angParam = isDirectionParam(partial1) ? hIdx1 : hIdx2;
        const auto posParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial2 : partial1);
        // clang-format off
        m_hessian[resIdx] = -m_gradProjDir[angParam].cross(wireDir).dot(line.gradient(posParam)) * Vector3::Unit(bending);
        // clang-format on
      }
      ACTS_VERBOSE("updateStrawResidual() - Hessian of the residual is "
                   << m_hessian[resIdx][bending]);
    }
  }
  return true;
}
bool StrawLineFitAuxiliaries::updateStrawAuxiliaries(const Line_t& line,
                                                     const Vector& wireDir) {
  const Vector& lineDir = line.direction();
  /// Between two calls the wire projection has not changed
  const double wireProject = lineDir.dot(wireDir);
  constexpr double s_tolerance = 1.e-12;

  if (false && std::abs(wireProject - m_wireProject) < s_tolerance) {
    ACTS_VERBOSE("Projection of the line matches the previous one."
                 << " Don't update the auxiliaries");
    return m_invProjDirLenSq > s_tolerance;
  }
  m_wireProject = wireProject;
  const double projDirLenSq = 1. - square(m_wireProject);
  /// The line is parallel to the wire
  if (projDirLenSq < s_tolerance) {
    ACTS_VERBOSE("Line & wire are parallel: " << toString(wireDir) << " vs. "
                                              << toString(lineDir));
    m_invProjDirLenSq = 0.;
    reset();
    return false;
  }
  m_invProjDirLenSq = 1. / projDirLenSq;
  m_invProjDirLen = std::sqrt(m_invProjDirLenSq);
  /// Project the segment line onto the wire plane and normalize
  m_projDir = (lineDir - m_wireProject * wireDir) * m_invProjDirLen;
  /// Loop over all configured parameters and calculate the partials
  /// of the wire projection
  for (auto partial : m_cfg.parsToUse) {
    /// Skip the parameters that are not directional
    if (!isDirectionParam(partial)) {
      continue;  // skip these parameters
    }
    const std::size_t pIdx = static_cast<std::size_t>(partial);
    const Vector& lGrad{line.gradient(static_cast<LineIndex>(partial))};
    m_projDirLenPartial[pIdx] = lGrad.dot(wireDir);
    // clang-format off
    m_gradProjDir[pIdx] = m_invProjDirLen * (lGrad - m_projDirLenPartial[pIdx] * wireDir) +
                           m_projDirLenPartial[pIdx] * m_wireProject * m_projDir * m_invProjDirLenSq;
    // clang-format on
    if (!m_cfg.useHessian) {
      continue;  // skip the Hessian calculation
    }
    for (auto partial1 : m_cfg.parsToUse) {
      /// Skip the parameters that are not directional or avoid double counting
      if (!isDirectionParam(partial1)) {
        continue;
      } else if (partial1 > partial) {
        break;
      }
      const auto param = static_cast<std::size_t>(partial);
      const auto param1 = static_cast<std::size_t>(partial1);
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

void StrawLineFitAuxiliaries::updateAlongTheStraw(const Line_t& line,
                                                  const Vector& hitMinSeg,
                                                  const Vector& wireDir) {
  m_residual[nonBending] =
      m_invProjDirLenSq * (hitMinSeg.dot(line.direction()) * m_wireProject -
                           hitMinSeg.dot(wireDir));
  ACTS_VERBOSE("updateAlongTheStraw() - Residual is "
               << m_residual[nonBending]);

  for (const auto partial : m_cfg.parsToUse) {
    if (partial == FitParIndex::t0) {
      continue;
    }
    const auto param = static_cast<std::size_t>(partial);
    const Vector& lGradient = line.gradient(static_cast<LineIndex>(partial));
    if (isDirectionParam(partial)) {
      // clang-format off
      m_gradient[param][nonBending] = m_invProjDirLenSq * (
                                       hitMinSeg.dot(lGradient) * m_wireProject +
                                       hitMinSeg.dot(line.direction()) * m_projDirLenPartial[param] +
                                       2. * m_residual[nonBending] * m_wireProject * m_projDirLenPartial[param]);
      // clang-format on

    } else if (isPositionParam(partial)) {
      // clang-format off
      m_gradient[param][nonBending] = -m_invProjDirLenSq * (lGradient.dot(line.direction()* m_wireProject - wireDir));
      // clang-format on
    }
    ACTS_VERBOSE("updateAlongTheStraw() - First derivative w.r.t "
                 << parName(partial) << " is "
                 << m_gradient[param][nonBending]);
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
      if (!(isDirectionParam(partial2) || isDirectionParam(partial1))) {
        continue;
      }
      const auto param = static_cast<std::size_t>(partial1);
      const auto param1 = static_cast<std::size_t>(partial2);
      const auto hessIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        auto calcMixedHessian = [this, &hitMinSeg, &line](const std::size_t p1,
                                                          const std::size_t p2) -> double {
          const double term1 = hitMinSeg.dot(line.gradient(static_cast<LineIndex>(p1))) * m_projDirLenPartial[p2];
          const double term2 = m_residual[nonBending] * m_projDirLenPartial[p1] * m_projDirLenPartial[p2];
          const double term3 = 2. * m_wireProject * m_projDirLenPartial[p1] * m_gradient[p2][nonBending];
          return (term1 + term2 + term3) * m_invProjDirLenSq;
        };
        const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial1), static_cast<LineIndex>(partial2));
        const double projHess = lHessian.dot(wireDir);
        m_hessian[hessIdx][nonBending] = (param != param1 ? calcMixedHessian(param, param1) + calcMixedHessian(param1, param)
                                                          : 2. * calcMixedHessian(param1, param)) +
                                        m_invProjDirLenSq * (hitMinSeg.dot(lHessian) * m_wireProject +
                                                             hitMinSeg.dot(line.direction()) * projHess +
                                                             2. * m_residual[nonBending] * m_wireProject * projHess);
        // clang-format on
      } else {
        /// Mixed positional and angular derivative
        const auto angParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial1 : partial2);
        const auto posParam = static_cast<LineIndex>(
            isDirectionParam(partial1) ? partial2 : partial1);
        const auto angIdx = static_cast<std::size_t>(angParam);
        // clang-format off
        m_hessian[hessIdx][nonBending] = m_invProjDirLenSq * (
                -line.gradient(posParam).dot(line.gradient(angParam)) * m_wireProject -
                 line.gradient(posParam).dot(line.direction()) * m_projDirLenPartial[angIdx] +
                 2. * m_gradient[static_cast<std::size_t>(posParam)][nonBending] * m_wireProject * m_projDirLenPartial[angIdx]);
        // clang-format on
      }
      ACTS_VERBOSE("updateAlongTheStraw() - Second derivative w.r.t. "
                   << parName(partial1) << ", " << parName(partial2) << " is. "
                   << m_hessian[hessIdx][nonBending]);
    }
  }
}

void StrawLineFitAuxiliaries::updateStripResidual(
    const Line_t& line, const Vector& normal, const Vector& sensorN,
    const Vector& sensorD, const Vector& stripPos, const bool isBending,
    const bool isNonBending) {
  using namespace Acts::UnitLiterals;

  const double normDot = normal.dot(line.direction());

  constexpr double tolerance = 1.e-12;
  if (std::abs(normDot) < tolerance) {
    reset();
    ACTS_WARNING("Segment line is embedded into the strip plane "
                 << toString(line.direction())
                 << ", normal: " << toString(normal));
    return;
  }
  const double planeOffSet = normal.dot(stripPos);
  ACTS_VERBOSE("Plane normal: "
               << toString(normal) << " |" << normal.norm()
               << "|, line direction: " << toString(line.direction())
               << ", angle: " << angle(normal, line.direction()));
  const double travelledDist =
      (planeOffSet - line.position().dot(normal)) / normDot;
  ACTS_VERBOSE("Need to travel " << travelledDist
                                 << " mm to reach the strip plane. ");

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
  /// Linear independent vectors
  /// Normal vector is indeed normal onto the plane spaned by these two vectors
  assert(Acts::abs(b1.dot(normal)) < tolerance);
  assert(Acts::abs(b2.dot(normal)) < tolerance);
  assert(isBending || isNonBending);
  const bool decompStereo =
      (m_cfg.calcAlongStrip || (isBending && isNonBending)) &&
      Acts::abs(b1DotB2) > tolerance;
  if (decompStereo) {
    /// A = lambda B1 + kappa B2
    ///                <A, B2> - <A,B1><B1,B2>
    ///   --> kappa = --------------------------
    ///                    1 - <B1,B2>^{2}
    ///                <A, B1> - <A,B2><B1,B2>
    ///   --> lambda = --------------------------
    ///                    1 - <B1,B2>^{2}
    const double invDist = 1. / (1. - square(b1DotB2));
    m_stereoTrf(bending, bending) = m_stereoTrf(nonBending, nonBending) =
        invDist;
    m_stereoTrf(bending, nonBending) = m_stereoTrf(nonBending, bending) =
        -b1DotB2 * invDist;
  }
  ACTS_VERBOSE("Meausres bending "
               << (isBending ? "yes" : "no")
               << ", measures non-bending:" << (isNonBending ? "yes" : "no"));
  ACTS_VERBOSE("strip orientation b1: "
               << toString(b1) << " |" << b1.norm() << "|"
               << ", b2: " << toString(b2) << " |" << b2.norm() << "|"
               << ", stereo angle: " << angle(b1, b2));

  auto assignResidual = [&](const Vector& calcDistance, Vector& residual) {
    residual[bending] =
        isBending || m_cfg.calcAlongStrip ? b1.dot(calcDistance) : 0.;
    residual[nonBending] =
        isNonBending || m_cfg.calcAlongStrip ? b2.dot(calcDistance) : 0.;
    if (decompStereo) {
      auto spatial = residual.block<2, 1>(0, 0);
      spatial = m_stereoTrf * spatial;
    }
    residual[time] = 0.;
    ACTS_VERBOSE("Distance: " << toString(calcDistance)
                              << ", <calc, n> =" << calcDistance.dot(normal)
                              << " -> residual: " << toString(residual)
                              << " --> closure test:"
                              << toString(residual[bending] * b1 +
                                          residual[nonBending] * b2));
  };
  /// Update the residual accordingly
  assignResidual(line.position() + travelledDist * line.direction() - stripPos,
                 m_residual);
  for (const auto partial : m_cfg.parsToUse) {
    ACTS_VERBOSE("stripResidual() - Calculate partial derivative w.r.t "
                 << parName(partial));
    switch (partial) {
      case FitParIndex::phi:
      case FitParIndex::theta: {
        // clang-format off
        const Vector& lGrad = line.gradient(static_cast<LineIndex>(partial));
        const double partialDist = -travelledDist / normDot * normal.dot(lGrad);
        const Vector gradCmp = travelledDist * lGrad +
                               partialDist * line.direction();
        // clang-format on
        assignResidual(gradCmp, m_gradient[static_cast<std::uint8_t>(partial)]);
        break;
      }
      case FitParIndex::y0:
      case FitParIndex::x0: {
        // clang-format off
        const Vector& lGrad = line.gradient(static_cast<LineIndex>(partial));
        const Vector gradCmp = lGrad - lGrad.dot(normal) / normDot * line.direction();
        assignResidual(gradCmp, m_gradient[static_cast<std::uint8_t>(partial)]);
        // clang-format on
        break;
      }
      default:
        /// Don't calculate anything for time
        break;
    }
  }

  if (!m_cfg.useHessian) {
    return;
  }
  for (const auto partial1 : m_cfg.parsToUse) {
    for (const auto partial2 : m_cfg.parsToUse) {
      /// Avoid double counting
      if (partial1 > partial2) {
        break;
      }
      /// At least one parameter needs to be a directional one
      if (!(isDirectionParam(partial1) || isDirectionParam(partial2))) {
        continue;
      }
      ACTS_VERBOSE(
          "stripResidual() - Calculate second partial derivative w.r.t "
          << parName(partial1) << ", " << parName(partial2));
      const auto param = static_cast<std::size_t>(partial1);
      const auto param1 = static_cast<std::size_t>(partial2);
      const auto resIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(partial1) && isDirectionParam(partial2)) {
        // clang-format off
        auto calcMixedTerms = [&line, &normal, &normDot, &b1, &b2, this](const std::size_t p1,
                                                                         const std::size_t p2) {

          return -normal.dot(line.gradient(static_cast<LineIndex>(p1))) / normDot *
                          (m_gradient[p2][bending] * b1 + m_gradient[p2][nonBending] * b2);
        };
        const Vector& lHessian = line.hessian(static_cast<LineIndex>(partial1),
                                              static_cast<LineIndex>(partial2));
        const Vector hessianCmp = travelledDist * (lHessian - normal.dot(lHessian) / normDot * line.direction()) +
                                  calcMixedTerms(param1, param) + calcMixedTerms(param, param1);
        assignResidual(hessianCmp, m_hessian[resIdx]);
        // clang-format on
      } else {
        // clang-format off
        const auto angParam = static_cast<LineIndex>(isDirectionParam(partial1) ? partial1 : partial2);
        const auto posParam = static_cast<LineIndex>(isDirectionParam(partial1) ? partial2 : partial1);
        const double lH = (normal.dot(line.gradient(posParam)) / normDot);

        const Vector hessianCmp = lH * ((normal.dot(line.gradient(angParam)) / normDot) * line.direction() -
                                        line.gradient(angParam));
        // clang-format on
        assignResidual(hessianCmp, m_hessian[resIdx]);
      }
    }
  }
}

}  // namespace Acts::Experimental::detail
