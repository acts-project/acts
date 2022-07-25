// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "Acts/Seeding/SeedFinderUtils.hpp"

namespace Acts {
LinCircle transformCoordinates(Acts::SpacePoint& sp,
                               Acts::SpacePoint& spM,
                               bool bottom) {
  // The computation inside this function is exactly identical to that in the
  // vectorized version of this function, except that it operates on a single
  // spacepoint. Please see the other version of this function for more
  // detailed comments.

  float cosPhiM = spM.x() / spM.radius();
  float sinPhiM = spM.y() / spM.radius();
  float deltaX = sp.x() - spM.x();
  float deltaY = sp.y() - spM.y();
  float deltaZ = sp.z() - spM.z();
  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);
  int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
  float cot_theta = deltaZ * iDeltaR * bottomFactor;
  LinCircle l;
  l.cotTheta = cot_theta;
  l.Zo = spM.z() - spM.radius() * cot_theta;
  l.iDeltaR = iDeltaR;
  l.U = x * iDeltaR2;
  l.V = y * iDeltaR2;
  l.Er = ((spM.varianceZ() + sp.varianceZ()) +
          (cot_theta * cot_theta) * (spM.varianceR() + sp.varianceR())) *
         iDeltaR2;
  return l;
}

void transformCoordinates(std::vector<Acts::SpacePoint*>& vec,
                          Acts::SpacePoint& spM, bool bottom,
                          std::vector<LinCircle>& linCircleVec) {

  float cosPhiM = spM.x() / spM.radius();
  float sinPhiM = spM.y() / spM.radius();
  for (auto sp : vec) {

    float deltaX = sp->x() - spM.x();
    float deltaY = sp->y() - spM.y();
    float deltaZ = sp->z() - spM.z();
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY * sinPhiM;
    float y = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ * iDeltaR * bottomFactor;
    // VERY frequent O(SP^3) access
    LinCircle l;
    l.cotTheta = cot_theta;
    // location on z-axis of this SP-duplet
    l.Zo = spM.z() - spM.radius() * cot_theta;
    l.iDeltaR = iDeltaR;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    l.U = x * iDeltaR2;
    l.V = y * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    l.Er = ((spM.varianceZ() + sp->varianceZ()) +
            (cot_theta * cot_theta) * (spM.varianceR() + sp->varianceR())) *
           iDeltaR2;

    l.x = x;
    l.y = y;
    l.z = sp->z();
    l.r = sp->radius();

    linCircleVec.push_back(l);
    sp->setCotTheta(cot_theta);

    sp->setDeltaR(std::sqrt((x * x) + (y * y) + (deltaZ * deltaZ)));
  }
  // sort the SP in order of cotTheta
  std::sort(vec.begin(), vec.end(),
            [](Acts::SpacePoint* a, Acts::SpacePoint* b) -> bool {
              return (a->cotTheta() < b->cotTheta());
            });
  std::sort(linCircleVec.begin(), linCircleVec.end(),
            [](const LinCircle& a, const LinCircle& b) -> bool {
              return (a.cotTheta < b.cotTheta);
            });
}

bool xyzCoordinateCheck(Acts::SeedfinderConfig m_config,
                        Acts::SpacePoint* sp,
                        const double* spacepointPosition,
                        const float toleranceParam,
                        double* outputCoordinates) {
  // check the compatibility of SPs coordinates in xyz assuming the
  // Bottom-Middle direction with the strip measurement details

  const float topHalfStripLength = m_config.getTopHalfStripLength(*sp);
  const float bottomHalfStripLength =
      m_config.getBottomHalfStripLength(*sp);
  const Acts::Vector3 topStripDirection =
      m_config.getTopStripDirection(*sp);
  const Acts::Vector3 bottomStripDirection =
      m_config.getBottomStripDirection(*sp);
  const Acts::Vector3 stripCenterDistance =
      m_config.getStripCenterDistance(*sp);

  // cross product between top strip vector and spacepointPosition
  double d1[3] = {
      (topHalfStripLength * topStripDirection[1]) * spacepointPosition[2] -
          (topHalfStripLength * topStripDirection[2]) * spacepointPosition[1],
      (topHalfStripLength * topStripDirection[2]) * spacepointPosition[0] -
          (topHalfStripLength * topStripDirection[0]) * spacepointPosition[2],
      (topHalfStripLength * topStripDirection[0]) * spacepointPosition[1] -
          (topHalfStripLength * topStripDirection[1]) * spacepointPosition[0]};

  // scalar product between bottom strip vector and d1
  double bd1 = (bottomHalfStripLength * bottomStripDirection[0]) * d1[0] +
               (bottomHalfStripLength * bottomStripDirection[1]) * d1[1] +
               (bottomHalfStripLength * bottomStripDirection[2]) * d1[2];

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  double s1 = (stripCenterDistance[0] * d1[0] + stripCenterDistance[1] * d1[1] +
               stripCenterDistance[2] * d1[2]);
  if (std::abs(s1) > std::abs(bd1) * toleranceParam)
    return false;

  // cross product between bottom strip vector and spacepointPosition
  double d0[3] = {(bottomHalfStripLength * bottomStripDirection[1]) *
                          spacepointPosition[2] -
                      (bottomHalfStripLength * bottomStripDirection[2]) *
                          spacepointPosition[1],
                  (bottomHalfStripLength * bottomStripDirection[2]) *
                          spacepointPosition[0] -
                      (bottomHalfStripLength * bottomStripDirection[0]) *
                          spacepointPosition[2],
                  (bottomHalfStripLength * bottomStripDirection[0]) *
                          spacepointPosition[1] -
                      (bottomHalfStripLength * bottomStripDirection[1]) *
                          spacepointPosition[0]};

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  double s0 = (stripCenterDistance[0] * d0[0] + stripCenterDistance[1] * d0[1] +
               stripCenterDistance[2] * d0[2]);
  if (std::abs(s0) > std::abs(bd1) * toleranceParam)
    return false;

  // if arive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Acts::Vector3 topStripCenterPosition =
      m_config.getTopStripCenterPosition(*sp);

  // spacepointPosition corected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates[0] = topStripCenterPosition[0] +
                         (topHalfStripLength * topStripDirection[0]) * s0;
  outputCoordinates[1] = topStripCenterPosition[1] +
                         (topHalfStripLength * topStripDirection[1]) * s0;
  outputCoordinates[2] = topStripCenterPosition[2] +
                         (topHalfStripLength * topStripDirection[2]) * s0;
  return true;
}
}  // namespace Acts
