// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/doublet.hpp"
#include "traccc/seeding/detail/lin_circle.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_type.hpp"

namespace traccc {

// helper functions used for both cpu and gpu
struct doublet_finding_helper {
  /// Check if two spacepoints form doublets
  ///
  /// @param sp1 is middle spacepoint
  /// @param sp2 is bottom or top spacepoint
  /// @param config is configuration parameter, with its fallbacks resolved
  ///        (@c seedfinder_config::resolve_fallbacks)
  /// @tparam otherSpType is whether it is for middle-bottom or middle-top
  /// doublet
  /// @tparam use_dphi_cut is whether to apply the azimuthal separation cut
  ///
  /// @return boolean value for compatibility
  ///
  template <details::spacepoint_type otherSpType, bool use_dphi_cut,
            typename T1, typename T2>
  static inline TRACCC_HOST_DEVICE bool isCompatible(
      const edm::spacepoint<T1>& sp1, const edm::spacepoint<T2>& sp2,
      const seedfinder_config& config);

  /// Check the azimuthal separation of two spacepoints against the maximum
  /// separation that a track with an impact parameter up to
  /// @c doubletDPhiD0Max can produce between their radii
  ///
  /// The bound is b = doubletDPhiConst + doubletDPhiSlope * |r2 - r1| +
  /// min(doubletDPhiCap, |asin(d0 / r2) - asin(d0 / r1)|), and the check
  /// |dphi| <= b is evaluated as sin(b) (p1 . p2) >= cos(b) |p1 x p2|, with
  /// sin(b) and cos(b) built from angle sums, so that no inverse
  /// trigonometric function is needed. The configuration terms are assumed
  /// to be non-negative.
  ///
  /// @return boolean value for compatibility
  ///
  static inline TRACCC_HOST_DEVICE bool isCompatibleDPhi(
      scalar x1, scalar y1, scalar r1, scalar x2, scalar y2, scalar r2,
      const seedfinder_config& config);

  /// Do the conformal transformation on doublet's coordinate
  ///
  /// @param sp1 is middle spacepoint
  /// @param sp2 is bottom or top spacepoint
  /// @tparam otherSpType is whether it is for middle-bottom or middle-top
  /// doublet
  ///
  /// @return lin_circle which contains the transformed coordinate information
  ///
  template <details::spacepoint_type otherSpType, typename T1, typename T2>
  static inline TRACCC_HOST_DEVICE lin_circle transform_coordinates(
      const edm::spacepoint<T1>& sp1, const edm::spacepoint<T2>& sp2);
};

template <details::spacepoint_type otherSpType, bool use_dphi_cut, typename T1,
          typename T2>
bool TRACCC_HOST_DEVICE doublet_finding_helper::isCompatible(
    const edm::spacepoint<T1>& sp1, const edm::spacepoint<T2>& sp2,
    const seedfinder_config& config) {
  static_assert(otherSpType == details::spacepoint_type::bottom ||
                otherSpType == details::spacepoint_type::top);

  const scalar r1 = sp1.radius();
  const scalar r2 = sp2.radius();

  scalar deltaR, cotTheta, zOrigin;
  if constexpr (otherSpType == details::spacepoint_type::bottom) {
    // check if R distance is too small, because bins are not R-sorted
    deltaR = r1 - r2;
    // actually cotTheta * deltaR to avoid division by 0 statements
    cotTheta = sp1.z() - sp2.z();
    // actually zOrigin * deltaR to avoid division by 0 statements
    zOrigin = sp1.z() * deltaR - r1 * cotTheta;
  } else {
    // check if R distance is too small, because bins are not R-sorted
    deltaR = r2 - r1;
    // actually cotTheta * deltaR to avoid division by 0 statements
    cotTheta = (sp2.z() - sp1.z());
    // actually zOrigin * deltaR to avoid division by 0 statements
    zOrigin = sp1.z() * deltaR - r1 * cotTheta;
  }

  scalar deltaRMin, deltaRMax;
  if constexpr (otherSpType == details::spacepoint_type::bottom) {
    deltaRMin = config.deltaRMinBottomSP;
    deltaRMax = config.deltaRMaxBottomSP;
  } else {
    deltaRMin = config.deltaRMinTopSP;
    deltaRMax = config.deltaRMaxTopSP;
  }

  if ((deltaR >= deltaRMax) || (deltaR <= deltaRMin) ||
      (math::fabs(cotTheta) >= config.cotThetaMax * deltaR) ||
      (zOrigin <= config.collisionRegionMin * deltaR) ||
      (zOrigin >= config.collisionRegionMax * deltaR) ||
      math::fabs(cotTheta) >= config.deltaZMax) {
    return false;
  }

  if constexpr (use_dphi_cut) {
    if (!isCompatibleDPhi(sp1.x(), sp1.y(), r1, sp2.x(), sp2.y(), r2, config)) {
      return false;
    }
  }

  if (!config.interactionPointCut) {
    return true;
  }

  /*
   * The following cut is capable of discriminating some doublets on the
   * basis that it is impossible to find a third spacepoint for the doublet
   * that will keep the resulting triplet inside the helix radius bound.
   * This explanation is enhanced with Geogebra commands with can be entered
   * into the application directly to provide a visual "proof" of why this
   * cut works.
   *
   * We will start by creating two spacepoints at arbitrary locations (they
   * can be moved as desired):
   *
   * ```
   * A = (2, 4)
   * B = (3, 12)
   * ```
   *
   * We will also define a radius $R$:
   *
   * ```
   * R = 10
   * ```
   *
   * Next, we consider the fact that two points and a radius define exactly
   * two circles through those points and with that radius. This makes sense
   * because three points (six degrees of freedom) precisely define a single
   * circle, and two points and a radius (five DoFs) define two circles. We
   * will construct those circles now, with the radius being the minimum
   * helix radius from the configuration.
   *
   * To find the midpoints of these circles, we will first find the
   * perpendicular bisector of points $A$ and $B$:
   *
   * ```
   * M = 0.5 * (A + B)
   * ```
   */
  scalar midX = 0.5f * (sp1.x() + sp2.x());
  scalar midY = 0.5f * (sp1.y() + sp2.y());

  /*
   * Then we will compute the slope of the perpendicular bisector:
   *
   * ```
   * s = (y(B) - y(A)) / (x(B) - x(A))
   * ```
   */
  scalar slope = (sp2.y() - sp1.y()) / (sp2.x() - sp1.x());

  /*
   * Next, we can simply place circle midpoints on our perpendicular
   * bisector, but we cannot simply use the radius $R$ as the distance from
   * the midpoint of $A$ and $B$, we have account for the length of the
   * sagitta:
   *
   * ```
   * dX = x(B) - x(A)
   * dY = y(B) - y(A)
   * dXY2 = dX * dX + dY * dY
   * q = sqrt(R * R - dXY2 / 4)
   * ```
   */
  scalar deltaX = sp2.x() - sp1.x();
  scalar deltaY = sp2.y() - sp1.y();
  scalar deltaXY2 = deltaX * deltaX + deltaY * deltaY;
  scalar sagittaLength = math::sqrt(
      config.minHelixRadius * config.minHelixRadius - deltaXY2 / 4.f);

  /*
   * We then compute the central angle between the points $A$, $B$, and the
   * midpoints of the circles we want to construct. Naively, this can be
   * done using the trigonomic functions:
   *
   * ```
   * theta = atan(1 / slope)
   * ```
   *
   * After which we can compute the delta between the midpoint of $A$ and
   * $B$ and the midpoints of our circles:
   *
   * ```
   * mdX = (R - q) * cos(theta)
   * mdY = (R - q) * sin(theta)
   * ```
   *
   * However, some trigonomy allows us to reduce this:
   *
   * ```
   * denom = sqrt((s * s + 1) / (s * s))
   * cosTheta = 1 / denom
   * sinTheta = -1 / (s * denom)
   * mdX = q * cosTheta
   * mdY = q * sinTheta
   * ```
   */
  scalar denom = math::sqrt((slope * slope + 1) / (slope * slope));
  scalar cosCentralAngle = 1.f / denom;
  scalar sinCentralAngle = -1.f / (slope * denom);

  scalar mpDeltaX = sagittaLength * cosCentralAngle;
  scalar mpDeltaY = sagittaLength * sinCentralAngle;

  /*
   * We now easily find the midpoints of the required circles:
   *
   * ```
   * V = (x(M) + mdX, y(M) + mdY)
   * W = (x(M) - mdX, y(M) - mdY)
   * ```
   */
  scalar mp1X = midX + mpDeltaX;
  scalar mp2X = midX - mpDeltaX;
  scalar mp1Y = midY + mpDeltaY;
  scalar mp2Y = midY - mpDeltaY;

  /*
   * Finally, we compute the radii of the circles. The crucial intuition
   * here is that if either of the newly constructed circles _completely_
   * contain a circle of radius $R$ around the origin, then we can never
   * find a third spacepoint to complete the triplet. Thus, we compute the
   * radii. We leave them squared in the C++ code to avoid an unnecessary
   * square root.
   */
  scalar mp1R2 = mp1X * mp1X + mp1Y * mp1Y;
  scalar mp2R2 = mp2X * mp2X + mp2Y * mp2Y;

  if (math::min(mp1R2, mp2R2) <= ((config.minHelixRadius - config.impactMax) *
                                  (config.minHelixRadius - config.impactMax))) {
    return false;
  }

  return true;
}

bool TRACCC_HOST_DEVICE doublet_finding_helper::isCompatibleDPhi(
    const scalar x1, const scalar y1, const scalar r1, const scalar x2,
    const scalar y2, const scalar r2, const seedfinder_config& config) {
  // Linear part of the bound
  const scalar a =
      config.doubletDPhiConst + config.doubletDPhiSlope * math::fabs(r2 - r1);
  if (a >= constant<scalar>::pi) {
    return true;
  }
  const scalar sinA = math::sin(a);
  const scalar cosA = math::cos(a);

  // Swing s = |asin(v2) - asin(v1)| of a track with impact parameter d0,
  // through sin(asin(v)) = v and cos(asin(v)) = sqrt(1 - v^2). It lies in
  // [0, pi/2].
  const scalar d0 = config.doubletDPhiD0Max;
  const scalar v1 = math::min(1.f, d0 / math::max(r1, 1.f));
  const scalar v2 = math::min(1.f, d0 / math::max(r2, 1.f));
  const scalar c1 = math::sqrt(1.f - v1 * v1);
  const scalar c2 = math::sqrt(1.f - v2 * v2);
  scalar sinM = math::fabs(v2 * c1 - c2 * v1);
  scalar cosM = c1 * c2 + v1 * v2;

  // min(cap, s): on [0, pi/2] the sine rises and the cosine falls
  if (config.doubletDPhiCap < 0.5f * constant<scalar>::pi) {
    sinM = math::min(sinM, math::sin(config.doubletDPhiCap));
    cosM = math::max(cosM, math::cos(config.doubletDPhiCap));
  }

  // b = a + min(cap, s) lies in [a, a + pi/2]; the comparison below needs
  // b <= pi, and b >= pi accepts every separation
  const scalar sinB = sinA * cosM + cosA * sinM;
  const scalar cosB = cosA * cosM - sinA * sinM;
  if (a >= 0.5f * constant<scalar>::pi && sinB <= 0.f) {
    return true;
  }

  // r1 r2 cos(dphi) and r1 r2 |sin(dphi)|
  const scalar dot = x1 * x2 + y1 * y2;
  const scalar cross = math::fabs(x1 * y2 - y1 * x2);
  return sinB * dot >= cosB * cross;
}

template <details::spacepoint_type otherSpType, typename T1, typename T2>
lin_circle TRACCC_HOST_DEVICE doublet_finding_helper::transform_coordinates(
    const edm::spacepoint<T1>& sp1, const edm::spacepoint<T2>& sp2) {
  static_assert(otherSpType == details::spacepoint_type::bottom ||
                otherSpType == details::spacepoint_type::top);

  const scalar& xM = sp1.x();
  const scalar& yM = sp1.y();
  const scalar& zM = sp1.z();
  const scalar& rM = sp1.radius();
  const scalar& varianceZM = sp1.z_variance();
  const scalar& varianceRM = sp1.radius_variance();
  scalar cosPhiM = xM / rM;
  scalar sinPhiM = yM / rM;

  scalar deltaX = sp2.x() - xM;
  scalar deltaY = sp2.y() - yM;
  scalar deltaZ = sp2.z() - zM;
  // calculate projection fraction of spM->sp vector pointing in same
  // direction as
  // vector origin->spM (x) and projection fraction of spM->sp vector pointing
  // orthogonal to origin->spM (y)
  scalar x = deltaX * cosPhiM + deltaY * sinPhiM;
  scalar y = deltaY * cosPhiM - deltaX * sinPhiM;
  // 1/(length of M -> SP)
  scalar iDeltaR2 =
      static_cast<scalar>(1.) / (deltaX * deltaX + deltaY * deltaY);
  scalar iDeltaR = std::sqrt(iDeltaR2);
  // cot_theta = (deltaZ/deltaR)
  scalar cot_theta = deltaZ * iDeltaR;
  if constexpr (otherSpType == details::spacepoint_type::bottom) {
    cot_theta = -cot_theta;
  }
  // VERY frequent (SP^3) access
  lin_circle l;
  l.m_cotTheta = cot_theta;
  // location on z-axis of this SP-duplet
  l.m_Zo = zM - rM * cot_theta;
  l.m_iDeltaR = iDeltaR;
  // transformation of circle equation (x,y) into linear equation (u,v)
  // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
  // is transformed into
  // 1 - 2x_0*u - 2y_0*v = 0
  // using the following m_U and m_V
  // (u = A + B*v); A and B are created later on
  l.m_U = x * iDeltaR2;
  l.m_V = y * iDeltaR2;
  // error term for sp-pair without correlation of middle space point
  l.m_Er = ((varianceZM + sp2.z_variance()) +
            (cot_theta * cot_theta) * (varianceRM + sp2.radius_variance())) *
           iDeltaR2;

  return l;
}

}  // namespace traccc
