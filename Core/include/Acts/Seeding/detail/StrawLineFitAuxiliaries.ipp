// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Seeding/detail/StrawLineFitAuxiliaries.hpp"

#include "Acts/Utilities/StringHelpers.hpp"

namespace Acts::detail {

template <StationSpacePoint Point_t>
void StrawLineFitAuxiliaries::updateStrawResidual(const Line_t& line,
                                                  const Point_t& strawMeas) {
  /** Fetch the hit position & direction */
  const auto& hitDir{strawMeas.sensorDirection()};
  if (!updateStrawAuxiliaries(line, hitDir)) {
    ACTS_WARNING("The line "
                 << toString(line.direction())
                 << " is parallel to the measurement: " << toString(hitDir));
    return;
  }

  /** Calculate the distance from the two reference points  */
  const Vector hitMinSeg = strawMeas.localPosition() - line.position();
  /** Distance from the segment line to the tube wire */
  const double lineDist = m_projDir.cross(hitDir).dot(hitMinSeg);
  const double resVal = (lineDist - strawMeas.driftRadius());
  m_residual = resVal * Vector::Unit(bending);

  /// If the tube is a twin-tube, the hit position is no longer arbitrary along
  /// the wire. Calculate the distance along the wire towards the point of
  /// closest approach.
#ifdef STONJEK
  if (hit.dimension() == 2) {
    m_residual[nonBending] =
        (hitMinSeg.dot(line.direction() ection()) * m_wireProject -
         hitMinSeg.dot(hitDir)) *
        m_invProjLenSq;
  }
#endif
  /** Calculate the first derivative of the residual */
  for (const auto param : m_cfg.parsToUse) {
    if (isDirectionParam(param)) {
      const double partialDist =
          m_gradProjDir[param].cross(hitDir).dot(hitMinSeg);
      ACTS_VERBOSE("Projection partial ("
                   << param << "): " << toString(m_gradProjDir[param])
                   << ", residual grad: " << partialDist);
      m_gradient[param] = partialDist * Vector::Unit(bending);
#ifdef STONJEK
      if (hit.dimension() == 2) {
        m_gradient[param][nonBending] =
            (hitMinSeg.dot(line.gradient[param]) * m_wireProject +
             hitMinSeg.dot(line.direction()) * m_gradProjDir[param]) *
                m_invProjLenSq +
            2. * resObj.residual[nonBending] *
                (m_wireProject * m_gradProjDir[param]) * m_invProjLenSq;
      }
#endif

    } else if (isPositionParam(param)) {
      const double partialDist =
          -m_projDir.cross(hitDir).dot(line.gradient(param));
      ACTS_VERBOSE("Parameter: " << param
                                 << ", partial residual: " << partialDist);
      m_gradient[param] = partialDist * Vector3::Unit(bending);
#ifdef STONJEK

      if (hit.dimension() == 2) {
        m_gradient[param][nonBending] =
            -(line.gradient[param].dot(line.direction()) * m_wireProject -
              line.gradient[param].dot(hitDir)) *
            m_invProjLenSq;
      }
#endif
    }
  }

  if (!m_cfg.useHessian) {
    return;
  }
  /** Loop to include the second order derivatvies */
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
      /// Second derivative w.r.t to the position parameters is zero
      if (!(isDirectionParam(param) || isDirectionParam(param1))) {
        continue;
      }
      const int resIdx = vecIdxFromSymMat<s_nLinePars>(param, param1);
      /// Pure angular derivatives of the residual
      if (isDirectionParam(param) && isDirectionParam(param1)) {
        const double partialSqDist =
            m_hessianProjDir[resIdx].cross(hitDir).dot(hitMinSeg);
        m_hessian[resIdx] = partialSqDist * Vector3::Unit(bending);
#ifdef STONJEK
        if (hit.dimension() == 2) {
          const double partialSqAlongWire =
              2. * resObj.residual[nonBending] * m_wireProject *
                  partSqLineProject * m_invProjLenSq +
              2. * resObj.residual[nonBending] * m_gradProjDir[param] *
                  m_gradProjDir[param1] * m_invProjLenSq +
              2. * m_gradient[param1][nonBending] * m_wireProject *
                  m_gradProjDir[param] * m_invProjLenSq +
              2. * m_gradient[param][nonBending] * m_wireProject *
                  m_gradProjDir[param1] * m_invProjLenSq +
              hitMinSeg.dot(line.hessian(param, param1)) * m_wireProject *
                  m_invProjLenSq +
              hitMinSeg.dot(line.direction()) * partSqLineProject *
                  m_invProjLenSq +
              hitMinSeg.dot(line.gradient[param]) * m_gradProjDir[param1] *
                  m_invProjLenSq +
              hitMinSeg.dot(line.gradient(param1)) * m_gradProjDir[param] *
                  m_invProjLenSq;
          resObj.hessian[resIdx][nonBending] = partialSqAlongWire;
        }
#endif
        continue;
      }
      /// Angular & Spatial mixed terms
      const auto angParam = isDirectionParam(param) ? param : param1;
      const auto posParam = isDirectionParam(param) ? param1 : param;

      const double partialSqDist =
          -m_gradProjDir[angParam].cross(hitDir).dot(line.gradient(posParam));
      m_hessian[resIdx] = partialSqDist * Vector3::Unit(bending);
#ifdef STONJEK
      if (hit.dimension() == 2) {
        const double partialSqAlongWire =
            -(line.gradient(param1).dot(line.gradient[param]) * m_wireProject +
              line.gradient(param1).dot(line.direction()) *
                  m_gradProjDir[param]) *
                m_invProjLenSq +
            2. * m_gradient[param1][nonBending] *
                (m_wireProject * m_gradProjDir[param]) * m_invProjLenSq;
        resObj.hessian[resIdx][nonBending] = partialSqAlongWire;
      }
#endif
    }
  }
}

template <StationSpacePoint Point_t>
void StrawLineFitAuxiliaries::updateStripResidual(const Line_t& line,
                                                  const Point_t& stripMeas) {
  const Vector& hitPos = stripMeas.localPosition();
  const Vector& normal = stripMeas.planeNormal();
  const double planeOffSet = normal.dot(hitPos);

  const double normDot = normal.dot(line.direction());

  constexpr double tolerance = 1.e-12;
  if (std::abs(normDot) < tolerance) {
    reset();
    ACTS_WARNING("Segment line is embedded into the strip plane "
                 << toString(line.direction())
                 << ", normal: " << toString(normal));
    return;
  }

  const double travelledDist =
      (planeOffSet - line.position().dot(normal)) / normDot;

  const double rePlaneAngle =
      stripMeas.sensorDirection().dot(stripMeas.sensorNormal());
  /// Linear independent vectors
      assert( Acts::abs(rePlaneAngle -1. ) > tolerance );
  /// Normal vector is indeed normal onto the plane spaned by these two vectors
  assert( Acts::abs(stripMeas.sensorDirection().dot(normal)) < tolerance);
  assert( Acts::abs(stripMeas.sensorNormal().dot(normal)) < tolerance);
  
  auto orthoStripPlane = [&](const Vector& v1, const Vector& v2) {
    return Acts::abs(rePlaneAngle) < tolerance
               ? v1
               : (v1 - rePlaneAngle * v2) / (1. - square(rePlaneAngle));
  };

  const Vector3 stripDir =
      orthoStripPlane(stripMeas.sensorDirection(), stripMeas.sensorNormal());
  const Vector3 toNextStrip =
      orthoStripPlane(stripMeas.sensorNormal(), stripMeas.sensorDirection());

  auto assignResidual = [&stripDir, &toNextStrip](const Vector& calcDistance,
                                                  Vector& residual) {
    residual[bending] = toNextStrip.dot(calcDistance);
    residual[nonBending] = stripDir.dot(calcDistance);
  };
  /// Update the residual accordingly
  assignResidual(line.position() + travelledDist * line.direction() - hitPos,
                 m_residual);
  for (unsigned fitPar : m_cfg.parsToUse) {
    switch (fitPar) {
      case FitParIndices::phi:
      case FitParIndices::theta: {
        const double partialDist =
            -travelledDist / normDot * normal.dot(line.gradient(fitPar));
        assignResidual(travelledDist * line.gradient(fitPar) +
                           partialDist * line.direction(),
                       m_gradient[fitPar]);
        break;
      }
      case FitParIndices::y0:
      case FitParIndices::x0: {
        assignResidual(
            line.gradient(fitPar) -
                line.gradient(fitPar).dot(normal) / normDot * line.direction(),
            m_gradient[fitPar]);
        break;
      }
      default: {
        break;
      }
    }
  }

  if (!m_cfg.useHessian) {
    return;
  }

  for (const auto param : m_cfg.parsToUse) {
    for (const auto param1 : m_cfg.parsToUse) {
      if (param1 > param) {
        break;
      }
      if (!isDirectionParam(param) || !isDirectionParam(param1)) {
        continue;
      }
      const int resIdx = vecIdxFromSymMat<s_nPars>(param, param1);
      if (isDirectionParam(param) && isDirectionParam(param1)) {
        assignResidual(
            travelledDist * line.hessian(param, param1) -
                ///
                travelledDist * (normal.dot(line.hessian(param, param1))) /
                    normDot * line.direction() -
                ///
                normal.dot(line.gradient(param1)) / normDot *
                    m_gradient[param] -
                ///
                normal.dot(line.gradient(param)) / normDot * m_gradient[param1],
            ///
            ///
            m_hessian[resIdx]);
        continue;
      }
      const auto angParam = isDirectionParam(param) ? param : param1;
      const auto posParam = isDirectionParam(param) ? param1 : param;

      const double gradientDisplace = normal.dot(line.gradient(param1));
      if (gradientDisplace < tolerance) {
        m_hessian[resIdx].setZero();
        continue;
      }
      assignResidual(
          gradientDisplace * (normal.dot(line.gradient(angParam)) /
                                  square(normDot) * line.direction() -
                              line.gradient(posParam) / normDot),
          m_hessian[resIdx]);
    }
  }
}

}  // namespace Acts::detail
