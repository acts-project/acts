// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
namespace Acts::detail {

StrawLineFitAuxiliaries::StrawLineFitAuxiliaries(
    const Config& cfg, std::unique_ptr<const Logger> logger)
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

      m_hessianProjDir[idx] =
          (lHessian - partSqLineProject * wireDir) * m_invProjDirLen +
          ///
          (m_projDirLenPartial[param1] * m_wireProject) * m_invProjDirLenSq *
              m_gradProjDir[param] +
          ///
          (m_projDirLenPartial[param] * m_wireProject) * m_invProjDirLenSq *
              m_gradProjDir[param1] +
          ///
          (partSqLineProject * m_wireProject) * m_invProjDirLenSq * m_projDir +
          ///
          (m_projDirLenPartial[param1] * m_projDirLenPartial[param]) *
              square(m_invProjDirLenSq) * m_projDir;
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
  ACTS_VERBOSE("Residual along the straw is " << m_residual[nonBending]);

  for (const auto param : m_cfg.parsToUse) {
    if (isDirectionParam(param)) {
      m_gradient[param][nonBending] =
          (hitMinSeg.dot(line.gradient(param)) * m_wireProject +
           hitMinSeg.dot(line.direction()) * m_projDirLenPartial[param] +
           2. * m_residual[nonBending] * m_wireProject *
               m_projDirLenPartial[param]) *
          m_invProjDirLenSq;
      ACTS_VERBOSE("First derivative w.r.t " << param << " (directional) is "
                                             << m_gradient[param][nonBending]);
    } else if (isPositionParam(param)) {
      m_gradient[param][nonBending] =
          -(line.gradient(param).dot(line.direction()) * m_wireProject -
            line.gradient(param).dot(wireDir)) *
          m_invProjDirLenSq;
      ACTS_VERBOSE("First derivative w.r.t " << param << " (positional) is "
                                             << m_gradient[param][nonBending]);
    }
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
      if (!((isDirectionParam(param) || isDirectionParam(param1)))) {
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
            calcMixedHessian(param, param1) + calcMixedHessian(param1, param) +
            (hitMinSeg.dot(line.hessian(param, param1)) * m_wireProject +
             hitMinSeg.dot(line.direction()) * projHess +
             2. * m_residual[nonBending] * m_wireProject * projHess) *
                m_invProjDirLenSq;
        continue;
      }
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
  }
}

}  // namespace Acts::detail
