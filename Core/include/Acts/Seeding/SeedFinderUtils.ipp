// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
template <typename external_spacepoint_t>
inline LinCircle transformCoordinates(
    const InternalSpacePoint<external_spacepoint_t>& sp,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom) {
  auto extractFunction =
      [](const InternalSpacePoint<external_spacepoint_t>& obj)
      -> std::array<float, 6> {
    std::array<float, 6> output{obj.x(),      obj.y(),         obj.z(),
                                obj.radius(), obj.varianceR(), obj.varianceZ()};
    return output;
  };

  return transformCoordinates<InternalSpacePoint<external_spacepoint_t>>(
      sp, spM, bottom, std::move(extractFunction));
}

template <typename external_spacepoint_t, typename callable_t>
inline LinCircle transformCoordinates(const external_spacepoint_t& sp,
                                      const external_spacepoint_t& spM,
                                      bool bottom,
                                      callable_t&& extractFunction) {
  // The computation inside this function is exactly identical to that in the
  // vectorized version of this function, except that it operates on a single
  // spacepoint. Please see the other version of this function for more
  // detailed comments.

  auto [xM, yM, zM, rM, varianceRM, varianceZM] = extractFunction(spM);
  auto [xSP, ySP, zSP, rSP, varianceRSP, varianceZSP] = extractFunction(sp);

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  float deltaX = xSP - xM;
  float deltaY = ySP - yM;
  float deltaZ = zSP - zM;
  float xNewFrame = deltaX * cosPhiM + deltaY * sinPhiM;
  float yNewFrame = deltaY * cosPhiM - deltaX * sinPhiM;
  float deltaR2 = (xNewFrame * xNewFrame + yNewFrame * yNewFrame);
  float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
  float iDeltaR = std::sqrt(iDeltaR2);
  int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
  float cotTheta = deltaZ * iDeltaR * bottomFactor;

  // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
  // circle into straight lines in the u/v plane the line equation can
  // be described in terms of aCoef and bCoef, where v = aCoef * u +
  // bCoef
  const float U = xNewFrame * iDeltaR2;
  const float V = yNewFrame * iDeltaR2;

  // error term for sp-pair without correlation of middle space point
  const float Er = ((varianceZM + varianceZSP) +
                    (cotTheta * cotTheta) * (varianceRM + varianceRSP)) *
                   iDeltaR2;

  sp.setDeltaR(std::sqrt(deltaR2 + (deltaZ * deltaZ)));

  return fillLineCircle({cotTheta, iDeltaR, Er, U, V, xNewFrame, yNewFrame});
}

template <typename external_spacepoint_t>
inline void transformCoordinates(
    Acts::SpacePointData& spacePointData,
    const std::vector<InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<LinCircle>& linCircleVec) {
  auto extractFunction =
      [](const InternalSpacePoint<external_spacepoint_t>& obj)
      -> std::array<float, 6> {
    std::array<float, 6> output{obj.x(),      obj.y(),         obj.z(),
                                obj.radius(), obj.varianceR(), obj.varianceZ()};
    return output;
  };

  transformCoordinates<InternalSpacePoint<external_spacepoint_t>>(
      spacePointData, vec, spM, bottom, linCircleVec,
      std::move(extractFunction));
}

template <typename external_spacepoint_t, typename callable_t>
inline void transformCoordinates(Acts::SpacePointData& spacePointData,
                                 const std::vector<external_spacepoint_t*>& vec,
                                 const external_spacepoint_t& spM, bool bottom,
                                 std::vector<LinCircle>& linCircleVec,
                                 callable_t&& extractFunction) {
  std::vector<std::size_t> indexes(vec.size());
  for (unsigned int i(0); i < indexes.size(); i++) {
    indexes[i] = i;
  }

  auto [xM, yM, zM, rM, varianceRM, varianceZM] = extractFunction(spM);

  // resize + operator[] is faster then reserve and push_back
  linCircleVec.resize(vec.size());

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;

  for (std::size_t idx(0); idx < vec.size(); ++idx) {
    auto& sp = vec[idx];
    auto [xSP, ySP, zSP, rSP, varianceRSP, varianceZSP] = extractFunction(*sp);

    float deltaX = xSP - xM;
    float deltaY = ySP - yM;
    float deltaZ = zSP - zM;
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float xNewFrame = deltaX * cosPhiM + deltaY * sinPhiM;
    float yNewFrame = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    float deltaR2 = (xNewFrame * xNewFrame + yNewFrame * yNewFrame);
    float iDeltaR2 = 1. / deltaR2;
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cotTheta = deltaZ * iDeltaR * bottomFactor;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    float U = xNewFrame * iDeltaR2;
    float V = yNewFrame * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    float Er = ((varianceZM + sp->varianceZ()) +
                (cotTheta * cotTheta) * (varianceRM + sp->varianceR())) *
               iDeltaR2;

    linCircleVec[idx] =
        fillLineCircle({cotTheta, iDeltaR, Er, U, V, xNewFrame, yNewFrame});
    spacePointData.setDeltaR(sp->index(),
                             std::sqrt(deltaR2 + (deltaZ * deltaZ)));
  }
}

inline LinCircle fillLineCircle(
    const std::array<float, 7>& lineCircleVariables) {
  auto [cotTheta, iDeltaR, Er, U, V, xNewFrame, yNewFrame] =
      lineCircleVariables;

  // VERY frequent (SP^3) access
  LinCircle l{};
  l.cotTheta = cotTheta;
  l.iDeltaR = iDeltaR;
  l.U = U;
  l.V = V;
  l.Er = Er;
  l.x = xNewFrame;
  l.y = yNewFrame;

  return l;
}

template <typename external_spacepoint_t>
inline bool xyzCoordinateCheck(
    Acts::SpacePointData& spacePointData,
    const Acts::SeedFinderConfig<external_spacepoint_t>& m_config,
    const Acts::InternalSpacePoint<external_spacepoint_t>& sp,
    const double* spacepointPosition, double* outputCoordinates) {
  // check the compatibility of SPs coordinates in xyz assuming the
  // Bottom-Middle direction with the strip measurement details
  bool hasValueStored = spacePointData.hasDynamicVariable();
  if (not hasValueStored) {
    return false;
  }

  std::size_t index = sp.index();

  const float& topHalfStripLength = spacePointData.getTopHalfStripLength(index);
  const float& bottomHalfStripLength =
      spacePointData.getBottomHalfStripLength(index);
  const Acts::Vector3& topStripDirection =
      spacePointData.getTopStripDirection(index);
  const Acts::Vector3& bottomStripDirection =
      spacePointData.getBottomStripDirection(index);
  const Acts::Vector3& stripCenterDistance =
      spacePointData.getStripCenterDistance(index);

  // prepare variables
  double xTopStripVector = topHalfStripLength * topStripDirection[0];
  double yTopStripVector = topHalfStripLength * topStripDirection[1];
  double zTopStripVector = topHalfStripLength * topStripDirection[2];
  double xBottomStripVector = bottomHalfStripLength * bottomStripDirection[0];
  double yBottomStripVector = bottomHalfStripLength * bottomStripDirection[1];
  double zBottomStripVector = bottomHalfStripLength * bottomStripDirection[2];

  // cross product between top strip vector and spacepointPosition
  double d1[3] = {yTopStripVector * spacepointPosition[2] -
                      zTopStripVector * spacepointPosition[1],
                  zTopStripVector * spacepointPosition[0] -
                      xTopStripVector * spacepointPosition[2],
                  xTopStripVector * spacepointPosition[1] -
                      yTopStripVector * spacepointPosition[0]};

  // scalar product between bottom strip vector and d1
  double bd1 = xBottomStripVector * d1[0] + yBottomStripVector * d1[1] +
               zBottomStripVector * d1[2];

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  double s1 = (stripCenterDistance[0] * d1[0] + stripCenterDistance[1] * d1[1] +
               stripCenterDistance[2] * d1[2]);
  if (std::abs(s1) > std::abs(bd1) * m_config.toleranceParam) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  double d0[3] = {yBottomStripVector * spacepointPosition[2] -
                      zBottomStripVector * spacepointPosition[1],
                  zBottomStripVector * spacepointPosition[0] -
                      xBottomStripVector * spacepointPosition[2],
                  xBottomStripVector * spacepointPosition[1] -
                      yBottomStripVector * spacepointPosition[0]};

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  double s0 = (stripCenterDistance[0] * d0[0] + stripCenterDistance[1] * d0[1] +
               stripCenterDistance[2] * d0[2]);
  if (std::abs(s0) > std::abs(bd1) * m_config.toleranceParam) {
    return false;
  }

  // if arive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Acts::Vector3& topStripCenterPosition =
      spacePointData.getTopStripCenterPosition(index);

  // spacepointPosition corected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates[0] = topStripCenterPosition[0] + xTopStripVector * s0;
  outputCoordinates[1] = topStripCenterPosition[1] + yTopStripVector * s0;
  outputCoordinates[2] = topStripCenterPosition[2] + zTopStripVector * s0;
  return true;
}
}  // namespace Acts
