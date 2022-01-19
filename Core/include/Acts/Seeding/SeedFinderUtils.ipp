// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
template <typename external_spacepoint_t>
LinCircle transformCoordinates(
    const InternalSpacePoint<external_spacepoint_t>& sp,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom) {
  // The computation inside this function is exactly identical to that in the
  // vectorized version of this function, except that it operates on a single
  // spacepoint. Please see the other version of this function for more
  // detailed comments.
  float xM = spM.x();
  float yM = spM.y();
  float zM = spM.z();
  float rM = spM.radius();
  float varianceZM = spM.varianceZ();
  float varianceRM = spM.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  float deltaX = sp.x() - xM;
  float deltaY = sp.y() - yM;
  float deltaZ = sp.z() - zM;
  float x = deltaX * cosPhiM + deltaY * sinPhiM;
  float y = deltaY * cosPhiM - deltaX * sinPhiM;
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);
  int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
  float cot_theta = deltaZ * iDeltaR * bottomFactor;
  LinCircle l;
  l.cotTheta = cot_theta;
  l.Zo = zM - rM * cot_theta;
  l.iDeltaR = iDeltaR;
  l.U = x * iDeltaR2;
  l.V = y * iDeltaR2;
  l.Er = ((varianceZM + sp.varianceZ()) +
          (cot_theta * cot_theta) * (varianceRM + sp.varianceR())) *
         iDeltaR2;
  return l;
}

template <typename external_spacepoint_t>
void transformCoordinates(
    std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    bool enableCutsForSortedSP, std::vector<LinCircle>& linCircleVec) {
  float xM = spM.x();
  float yM = spM.y();
  float zM = spM.z();
  float rM = spM.radius();
  float varianceZM = spM.varianceZ();
  float varianceRM = spM.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  for (auto sp : vec) {
    float deltaX = sp->x() - xM;
    float deltaY = sp->y() - yM;
    float deltaZ = sp->z() - zM;
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
    // VERY frequent (SP^3) access
    LinCircle l;
    l.cotTheta = cot_theta;
    // location on z-axis of this SP-duplet
    l.Zo = zM - rM * cot_theta;
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
    l.Er = ((varianceZM + sp->varianceZ()) +
            (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
           iDeltaR2;
    linCircleVec.push_back(l);
  }
  // sort the SP in order of cotTheta
  if (enableCutsForSortedSP) {
    std::sort(linCircleVec.begin(), linCircleVec.end(),
              [](const LinCircle& a, const LinCircle& b) -> bool {
                return (a.cotTheta < b.cotTheta);
              });
  }
}
}  // namespace Acts
