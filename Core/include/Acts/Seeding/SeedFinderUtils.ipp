// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

namespace SeedFinderUtils {
template <typename callable_t, typename external_spacepoint_t>
using isSignatureCompatible = decltype(
    std::declval<callable_t&>() = std::declval<external_spacepoint_t>());
}

template <typename external_spacepoint_t>
LinCircle transformCoordinates(
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
      sp, spM, bottom, extractFunction);
}

template <typename external_spacepoint_t, typename callable_t>
LinCircle transformCoordinates(const external_spacepoint_t& sp,
                               const external_spacepoint_t& spM, bool bottom,
                               callable_t&& extractFunction) {
  using member_ptr_type =
      std::array<float, 6> (*)(const external_spacepoint_t&);

  static_assert(
      Concepts::is_detected<SeedFinderUtils::isSignatureCompatible,
                            member_ptr_type, decltype(extractFunction)>::value,
      "Callable given does not correspond exactly to required call signature "
      "-> std::array<float, 6>(const external_spacepoint_t&)");

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
  l.Er = ((varianceZM + varianceZSP) +
          (cot_theta * cot_theta) * (varianceRM + varianceRSP)) *
         iDeltaR2;
  return l;
}

template <typename external_spacepoint_t>
void transformCoordinates(
    const std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
    const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    bool enableCutsForSortedSP, std::vector<LinCircle>& linCircleVec) {
  auto extractFunction =
      [](const InternalSpacePoint<external_spacepoint_t>& obj)
      -> std::array<float, 6> {
    std::array<float, 6> output{obj.x(),      obj.y(),         obj.z(),
                                obj.radius(), obj.varianceR(), obj.varianceZ()};
    return output;
  };

  return transformCoordinates<InternalSpacePoint<external_spacepoint_t>>(
      vec, spM, bottom, enableCutsForSortedSP, linCircleVec, extractFunction);
}

template <typename external_spacepoint_t, typename callable_t>
void transformCoordinates(const std::vector<const external_spacepoint_t*>& vec,
                          const external_spacepoint_t& spM, bool bottom,
                          bool enableCutsForSortedSP,
                          std::vector<LinCircle>& linCircleVec,
                          callable_t&& extractFunction) {
  using member_ptr_type =
      std::array<float, 6> (*)(const external_spacepoint_t&);

  static_assert(
      Concepts::is_detected<SeedFinderUtils::isSignatureCompatible,
                            member_ptr_type, decltype(extractFunction)>::value,
      "Callable given does not correspond exactly to required call signature "
      "-> std::array<float, 6>(const external_spacepoint_t&)");

  auto [xM, yM, zM, rM, varianceRM, varianceZM] = extractFunction(spM);

  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  for (auto sp : vec) {
    auto [xSP, ySP, zSP, rSP, varianceRSP, varianceZSP] = extractFunction(*sp);

    float deltaX = xSP - xM;
    float deltaY = ySP - yM;
    float deltaZ = zSP - zM;
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
