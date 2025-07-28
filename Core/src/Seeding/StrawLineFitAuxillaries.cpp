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
namespace {

using Vector = StrawLineFitAuxiliaries::Vector;

std::string parName(const std::size_t idx) {
  using ParIdx = StrawLineFitAuxiliaries::FitParIndices;
  switch (idx) {
    using enum ParIdx;
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
  }
  return "unknown";
}

double angle(const Vector& v1, const Vector& v2) {
  using namespace Acts::UnitLiterals;
  return std::acos(std::clamp(v1.dot(v2), -1., 1.)) / 1_degree;
}
}  // namespace

StrawLineFitAuxiliaries::StrawLineFitAuxiliaries(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{cfg}, m_logger{std::move(logger)} {
  std::ranges::sort(m_cfg.parsToUse);
}

using Vector = StrawLineFitAuxiliaries::Vector;
const Vector& StrawLineFitAuxiliaries::residual() const {
  return m_residual;
}
const Vector& StrawLineFitAuxiliaries::gradient(const std::size_t par) const {
  assert(par < m_gradient.size());
  return m_gradient[par];
}
const Vector& StrawLineFitAuxiliaries::hessian(const std::size_t param,
                                               const std::size_t param1) const {
  const std::size_t idx = vecIdxFromSymMat<s_nPars>(param, param1);
  assert(idx < m_hesian.size());
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
  /** Fetch the hit position & direction */
  if (!updateStrawAuxiliaries(line, wireDir)) {
    ACTS_WARNING("updateStrawResidual() - The line "
                 << toString(line.direction())
                 << " is parallel to the measurement: " << toString(wireDir));
    return false;
  }
  /** Distance from the segment line to the tube wire */
  const double lineDist = m_projDir.cross(wireDir).dot(hitMinSeg);
  const double resVal = (lineDist - driftRadius);
  m_residual = resVal * Vector::Unit(bending);

  /** Calculate the first derivative of the residual */
  for (const auto param : m_cfg.parsToUse) {
    if (isDirectionParam(param)) {
      const double partialDist =
          m_gradProjDir[param].cross(wireDir).dot(hitMinSeg);
      ACTS_VERBOSE("updateStrawResidual() - Partial "
                   << parName(param) << ": " << toString(m_gradProjDir[param])
                   << ", residual grad: " << partialDist);
      m_gradient[param] = partialDist * Vector::Unit(bending);

    } else if (isPositionParam(param)) {
      const double partialDist =
          -m_projDir.cross(wireDir).dot(line.gradient(param));
      ACTS_VERBOSE("updateStrawResidual() - Partial " << parName(param) << ": "
                                                      << partialDist);
      m_gradient[param] = partialDist * Vector3::Unit(bending);
    }
  }

  if (!m_cfg.useHessian) {
    ACTS_VERBOSE("updateStrawResidual() - Skip Hessian calculation");
    return true;
  }
  /** Loop to include the second order derivatvies */
  for (const auto param : m_cfg.parsToUse) {
    if (param == FitParIndices::t0) {
      continue;
    }
    for (const auto param1 : m_cfg.parsToUse) {
      if (param1 == FitParIndices::t0) {
        continue;
      } else if (param1 > param) {
        break;
      }
      ACTS_VERBOSE("updateStrawResidual() - Calculate Hessian for parameters "
                   << parName(param) << ", " << parName(param1) << ".");
      /// Second derivative w.r.t to the position parameters is zero
      if (!(isDirectionParam(param) || isDirectionParam(param1))) {
        continue;
      }
      const int resIdx = vecIdxFromSymMat<s_nLinePars>(param, param1);
      /// Pure angular derivatives of the residual
      if (isDirectionParam(param) && isDirectionParam(param1)) {
        const double partialSqDist =
            m_hessianProjDir[resIdx].cross(wireDir).dot(hitMinSeg);
        m_hessian[resIdx] = partialSqDist * Vector3::Unit(bending);
      } else {
        /// Angular & Spatial mixed terms
        const auto angParam = isDirectionParam(param) ? param : param1;
        const auto posParam = isDirectionParam(param) ? param1 : param;
        m_hessian[resIdx] = -m_gradProjDir[angParam].cross(wireDir).dot(
                                line.gradient(posParam)) *
                            Vector3::Unit(bending);
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
    ACTS_VERBOSE(
        "Projection of the line matches the previous one. Don't update the "
        "auxiliaries");
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
  for (auto param : m_cfg.parsToUse) {
    /// Skip the parameters that are not directional
    if (!isDirectionParam(param)) {
      continue;  // skip these parameters
    }
    m_projDirLenPartial[param] = line.gradient(param).dot(wireDir);
    m_gradProjDir[param] =
        m_invProjDirLen *
            (line.gradient(param) - m_projDirLenPartial[param] * wireDir) +
        ///
        m_projDirLenPartial[param] * m_wireProject * m_projDir *
            m_invProjDirLenSq;

    if (!m_cfg.useHessian) {
      continue;  // skip the Hessian calculation
    }
    for (auto param1 : m_cfg.parsToUse) {
      /// Skip the parameters that are not directional or avoid double counting
      if (!isDirectionParam(param1)) {
        continue;
      } else if (param1 > param) {
        break;
      }
      const std::size_t idx = vecIdxFromSymMat<s_nLinePars>(param, param1);
      const Vector& lHessian = line.hessian(param, param1);
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
      m_hessianProjDir[idx] =
          (lHessian - partSqLineProject * wireDir) * m_invProjDirLen +
          ///
          m_invProjDirLenSq * (calcMixedHessian() +
                               (partSqLineProject * m_wireProject) * m_projDir);
      ACTS_VERBOSE("updateStrawAuxiliaries() - Hessian w.r.t. "
                   << parName(param1) << ", " << parName(param) << " is "
                   << toString(m_hessianProjDir[idx]));
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

  for (const auto param : m_cfg.parsToUse) {
    if (isDirectionParam(param)) {
      m_gradient[param][nonBending] =
          (hitMinSeg.dot(line.gradient(param)) * m_wireProject +
           hitMinSeg.dot(line.direction()) * m_projDirLenPartial[param] +
           2. * m_residual[nonBending] * m_wireProject *
               m_projDirLenPartial[param]) *
          m_invProjDirLenSq;

    } else if (isPositionParam(param)) {
      m_gradient[param][nonBending] =
          -(line.gradient(param).dot(line.direction()) * m_wireProject -
            line.gradient(param).dot(wireDir)) *
          m_invProjDirLenSq;
    }
    ACTS_VERBOSE("updateAlongTheStraw() - First derivative w.r.t "
                 << parName(param) << " is " << m_gradient[param][nonBending]);
  }
  if (!m_cfg.useHessian) {
    return;
  }

  for (auto param : m_cfg.parsToUse) {
    if (param == FitParIndices::t0) {
      continue;
    }
    for (auto param1 : m_cfg.parsToUse) {
      if (param1 == FitParIndices::t0) {
        continue;
      } else if (param1 > param) {
        break;
      }
      /// second derivative w.r.t intercepts is always 0
      if (!(isDirectionParam(param) || isDirectionParam(param1))) {
        continue;
      }
      const std::size_t hessIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(param) && isDirectionParam(param1)) {
        auto calcMixedHessian = [this, &hitMinSeg, &line](
                                    const std::size_t p1,
                                    const std::size_t p2) -> double {
          const double term1 =
              hitMinSeg.dot(line.gradient(p1)) * m_projDirLenPartial[p2];
          const double term2 = m_residual[nonBending] *
                               m_projDirLenPartial[p1] *
                               m_projDirLenPartial[p2];
          const double term3 = 2. * m_wireProject * m_projDirLenPartial[p1] *
                               m_gradient[p2][nonBending];
          return (term1 + term2 + term3) * m_invProjDirLenSq;
        };

        const double projHess = line.hessian(param1, param).dot(wireDir);
        m_hessian[hessIdx][nonBending] =
            (param != param1 ? calcMixedHessian(param, param1) +
                                   calcMixedHessian(param1, param)
                             : 2. * calcMixedHessian(param1, param)) +
            (hitMinSeg.dot(line.hessian(param, param1)) * m_wireProject +
             hitMinSeg.dot(line.direction()) * projHess +
             2. * m_residual[nonBending] * m_wireProject * projHess) *
                m_invProjDirLenSq;
      } else {
        /// Mixed positional and angular derivative
        const auto angParam = isDirectionParam(param) ? param : param1;
        const auto posParam = isDirectionParam(param) ? param1 : param;
        m_hessian[hessIdx][nonBending] =
            (-line.gradient(posParam).dot(line.gradient(angParam)) *
                 m_wireProject -
             ///
             line.gradient(posParam).dot(line.direction()) *
                 m_projDirLenPartial[angParam] +
             ///
             2. * m_gradient[posParam][nonBending] * m_wireProject *
                 m_projDirLenPartial[angParam]) *
            m_invProjDirLenSq;
      }
      ACTS_VERBOSE("updateAlongTheStraw() - Second derivative w.r.t. "
                   << parName(param) << ", " << parName(param1) << " is. "
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
  assert(measuresBending || measuresNonBending);
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
      residual.block<2, 1>(0, 0) = m_stereoTrf * residual.block<2, 1>(0, 0);
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
  for (const std::size_t fitPar : m_cfg.parsToUse) {
    ACTS_VERBOSE("stripResidual() - Calculate partial derivative w.r.t "
                 << parName(fitPar));
    switch (fitPar) {
      case FitParIndices::phi:
      case FitParIndices::theta: {
        const double partialDist =
            -travelledDist / normDot * normal.dot(line.gradient(fitPar));
        const Vector gradCmp = travelledDist * line.gradient(fitPar) +
                               partialDist * line.direction();
        assignResidual(gradCmp, m_gradient[fitPar]);
        break;
      }
      case FitParIndices::y0:
      case FitParIndices::x0: {
        const Vector gradCmp =
            line.gradient(fitPar) -
            line.gradient(fitPar).dot(normal) / normDot * line.direction();
        assignResidual(gradCmp, m_gradient[fitPar]);
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
  for (const auto param : m_cfg.parsToUse) {
    for (const auto param1 : m_cfg.parsToUse) {
      /// Avoid double counting
      if (param1 > param) {
        break;
      }
      /// At least one parameter needs to be a directional one
      if (!(isDirectionParam(param) || isDirectionParam(param1))) {
        continue;
      }
      ACTS_VERBOSE(
          "stripResidual() - Calculate second partial derivative w.r.t "
          << parName(param) << ", " << parName(param1));

      const int resIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(param) && isDirectionParam(param1)) {
        auto calcMixedTerms = [&line, &normal, &normDot, &b1, &b2, this](
                                  const std::size_t p1, const std::size_t p2) {
          return -normal.dot(line.gradient(p1)) / normDot *
                 (m_gradient[p2][bending] * b1 +
                  m_gradient[p2][nonBending] * b2);
        };
        const Vector hessianCmp =
            travelledDist * (line.hessian(param, param1) -
                             ///
                             normal.dot(line.hessian(param, param1)) / normDot *
                                 line.direction()) +
            ///
            calcMixedTerms(param1, param) + calcMixedTerms(param, param1);
        assignResidual(hessianCmp, m_hessian[resIdx]);
      } else {
        const auto angParam = isDirectionParam(param) ? param : param1;
        const auto posParam = isDirectionParam(param) ? param1 : param;
        const double lH = (normal.dot(line.gradient(posParam)) / normDot);
        const Vector hessianCmp =
            lH * ((normal.dot(line.gradient(angParam)) / normDot) *
                      line.direction() -
                  line.gradient(angParam));
        assignResidual(hessianCmp, m_hessian[resIdx]);
      }
    }
  }
}

}  // namespace Acts::Experimental::detail
