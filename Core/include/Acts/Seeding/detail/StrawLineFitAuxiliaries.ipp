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

namespace Acts::detail {


template <StationSpacePoint sp_t>
void StrawLineFitAuxiliaries::updateStrawResidual(const Line_t& line,
                                                  const sp_t& strawMeas) {
  /** Fetch the hit position & direction */
  const auto& hitDir{strawMeas.sensorDirection()};
  if (!updateStrawAuxillaries(line, hitDir)) {
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
        (hitMinSeg.dot(line.direction()) * m_wireProject -
         hitMinSeg.dot(hitDir)) *
        m_invProjLenSq;
  }
  #endif
  /** Calculate the first derivative of the residual */
  for (const auto param : m_cfg.parsToUse) {
    if (isDirectionParam(param)) {
     const double partialDist = m_gradProjDir[param].cross(hitDir).dot(hitMinSeg);
     m_gradient[param] = partialDist * Vector::Unit(bending);
#ifdef STONJEK
            if (hit.dimension() == 2) {
          m_gradient[param][nonBending] =
              (hitMinSeg.dot(line.gradient[param]) * m_wireProject +
               hitMinSeg.dot(line.dir) * m_gradProjDir[param]) *
                  m_invProjLenSq +
              2. * resObj.residual[nonBending] *
                  (m_wireProject * m_gradProjDir[param]) *
                  m_invProjLenSq;
        }
 #endif

    } else if(isPositionParam(param)) {

         const double partialDist =
            -m_projDir.cross(hitDir).dot(line.gradient(param));
        m_gradient[param] =
            partialDist * Vector3::Unit(bending);
    #ifdef STONJEK

            if (hit.dimension() == 2) {
          m_gradient[param][nonBending] =
              -(line.gradient[param].dot(line.dir) * m_wireProject -
                line.gradient[param].dot(hitDir)) *
              m_invProjLenSq;
        }
#endif   
    }
}

  if (!m_cfg.useHessian) {
    return;
  }
 #ifdef STONJEK
  /** Loop to include the second order derivatvies */
  for (int param = toInt(ParamDefs::phi); param >= 0; --param) {
    if (!resObj.evalPhiPars &&
        (param == toInt(ParamDefs::x0) || param == toInt(ParamDefs::phi))) {
      continue;
    }
    for (int param1 = param; param1 >= 0; --param1) {
      if (!resObj.evalPhiPars &&
          (param1 == toInt(ParamDefs::x0) || param1 == toInt(ParamDefs::phi))) {
        continue;
      }
      const int lineIdx = vecIdxFromSymMat<nLinePars>(param, param1);
      const int resIdx =
          vecIdxFromSymMat<toInt(ParamDefs::nPars)>(param, param1);
      /// Pure angular derivatives of the residual
      if ((param == toInt(ParamDefs::theta) ||
           param == toInt(ParamDefs::phi)) &&
          (param1 == toInt(ParamDefs::theta) ||
           param1 == toInt(ParamDefs::phi))) {
        const double partSqLineProject = line.hessian[lineIdx].dot(hitDir);
        const Vector3 projDirPartSq =
            (line.hessian[lineIdx] - partSqLineProject * hitDir) *
                m_invProjLen +
            (m_gradProjDir[param1] * m_wireProject) *
                m_invProjLenSq * m_gradProjDir[param] +
            (m_gradProjDir[param] * m_wireProject) *
                m_invProjLenSq * m_gradProjDir[param1] +
            (partSqLineProject * m_wireProject) *
                m_invProjLenSq * m_projDir +
            (m_gradProjDir[param1] *
             m_gradProjDir[param]) *
                std::pow(m_invProjLenSq, 2) * m_projDir;

        const double partialSqDist = projDirPartSq.cross(hitDir).dot(hitMinSeg);
        resObj.hessian[resIdx] =
            partialSqDist * Vector3::Unit(bending);
        if (hit.dimension() == 2) {
          const double partialSqAlongWire =
              2. * resObj.residual[nonBending] *
                  m_wireProject * partSqLineProject *
                  m_invProjLenSq +
              2. * resObj.residual[nonBending] *
                  m_gradProjDir[param] *
                  m_gradProjDir[param1] * m_invProjLenSq +
              2. * m_gradient[param1][nonBending] *
                  m_wireProject * m_gradProjDir[param] *
                  m_invProjLenSq +
              2. * m_gradient[param][nonBending] *
                  m_wireProject * m_gradProjDir[param1] *
                  m_invProjLenSq +
              hitMinSeg.dot(line.hessian[lineIdx]) * m_wireProject *
                  m_invProjLenSq +
              hitMinSeg.dot(line.dir) * partSqLineProject *
                  m_invProjLenSq +
              hitMinSeg.dot(line.gradient[param]) *
                  m_gradProjDir[param1] * m_invProjLenSq +
              hitMinSeg.dot(line.gradient[param1]) *
                  m_gradProjDir[param] * m_invProjLenSq;
          resObj.hessian[resIdx][nonBending] = partialSqAlongWire;
        }
      }
      /// Angular & Spatial mixed terms
      else if (param == toInt(ParamDefs::theta) ||
               param == toInt(ParamDefs::phi)) {
        const double partialSqDist =
            -m_gradProjDir[param].cross(hitDir).dot(line.gradient[param1]);
        resObj.hessian[resIdx] =
            partialSqDist * Vector3::Unit(bending);
        if (hit.dimension() == 2) {
          const double partialSqAlongWire =
              -(line.gradient[param1].dot(line.gradient[param]) *
                    m_wireProject +
                line.gradient[param1].dot(line.dir) *
                    m_gradProjDir[param]) *
                  m_invProjLenSq +
              2. * m_gradient[param1][nonBending] *
                  (m_wireProject * m_gradProjDir[param]) *
                  m_invProjLenSq;
          resObj.hessian[resIdx][nonBending] = partialSqAlongWire;
        }
      }
    }
  }
#endif
}

}  // namespace Acts::detail
