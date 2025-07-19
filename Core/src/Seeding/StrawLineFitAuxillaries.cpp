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
    : m_cfg{cfg}, m_logger{std::move(logger)} {}

using Vector = StrawLineFitAuxiliaries::Vector;
const Vector& StrawLineFitAuxiliaries::residual() const {
  return m_residual;
}
const Vector& StrawLineFitAuxiliaries::gradient(const unsigned par) const {
  assert(par < m_gradient.size());
  return m_gradient[par];
}
const Vector& StrawLineFitAuxiliaries::hessian(const unsigned param,
                                               const unsigned param1) const {
  const unsigned idx = vecIdxFromSymMat<s_nPars>(param, param1);
  assert(idx < m_hesian.size());
  return m_hessian[idx];
}
bool StrawLineFitAuxiliaries::updateStrawAuxiliaries(const Line_t& line,
                                                     const Vector& wireDir) {
  const Vector& lineDir = line.direction();
  const double wireProject = lineDir.dot(wireDir);
  /// Between two calls the wire projection has not changed
  if (std::abs(wireProject - m_wireProject) < s_tolerance) {
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
    m_residual.setZero();
    for (Vector3& grad : m_gradient) {
      grad.setZero();
    }
    if (m_cfg.useHessian) {
      for (Vector3& hess : m_hessian) {
        hess.setZero();
      }
    }
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
      if (!isDirectionParam(param1) || param1 > param) {
        continue;
      }
      const unsigned idx = vecIdxFromSymMat<s_nLinePars>(param, param1);
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
}  // namespace Acts::detail
