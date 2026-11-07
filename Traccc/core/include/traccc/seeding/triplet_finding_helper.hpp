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
#include "traccc/seeding/detail/strip_spacepoint_calibration.hpp"
#include "traccc/seeding/detail/triplet.hpp"

namespace traccc {

/// Quantities of a middle-bottom doublet needed by the strip aware triplet
/// finding, computed once per doublet
struct strip_triplet_context {
  scalar rM, varianceZM, varianceRM;
  // bottom doublet quantities
  scalar cotThetaB0, erB, iDeltaRB, Ub0, Vb0, xB0, yB0;
  scalar sinTheta, cosTheta;
  // x and y terms for the rotation from UV to XY plane
  scalar rotationTermsUVtoXY[2];
  scalar scatteringInRegion2, sigmaSquaredPtDependent;

  template <typename T>
  TRACCC_HOST_DEVICE strip_triplet_context(const edm::spacepoint<T>& spM,
                                           const lin_circle& lb,
                                           const seedfinder_config& config,
                                           const scalar iSinTheta2,
                                           const scalar scatteringInRegion2_)
      : rM(spM.radius()),
        varianceZM(spM.z_variance()),
        varianceRM(spM.radius_variance()),
        cotThetaB0(lb.cotTheta()),
        erB(lb.Er()),
        iDeltaRB(lb.iDeltaR()),
        Ub0(lb.U()),
        Vb0(lb.V()),
        xB0(lb.U() / (lb.iDeltaR() * lb.iDeltaR())),
        yB0(lb.V() / (lb.iDeltaR() * lb.iDeltaR())),
        sinTheta(1.f / std::sqrt(iSinTheta2)),
        cosTheta(lb.cotTheta() * sinTheta),
        rotationTermsUVtoXY{spM.x() / spM.radius() * sinTheta,
                            spM.y() / spM.radius() * sinTheta},
        scatteringInRegion2(scatteringInRegion2_),
        // 4 sigma^2 (highland / pTPerHelixRadius)^2 / sin^2(theta)
        sigmaSquaredPtDependent(iSinTheta2 * config.pT2perRadius * 4.f *
                                config.sigmaScattering *
                                config.sigmaScattering) {}
};

// helper function used for both cpu and gpu
struct triplet_finding_helper {
  /// Check if two doublets with common middle spacepoint can form a triplet,
  /// moving the strip spacepoints along their strips to the crossing of the
  /// triplet direction
  ///
  /// @param ctx is the middle-bottom doublet context
  /// @param calM is the calibration of the middle spacepoint
  /// @param calB is the calibration of the bottom spacepoint
  /// @param calT is the calibration of the top spacepoint
  /// @param lt is transformed coordinate of middle-top doublet
  /// @param config is configuration parameter
  /// @param curvature is curvature of triplet
  /// @param impact_parameter is impact parameter of triplet
  ///
  /// @return boolean value for compatibility
  static inline TRACCC_HOST_DEVICE bool isCompatibleStrip(
      const strip_triplet_context& ctx,
      const strip_spacepoint_calibration& calM,
      const strip_spacepoint_calibration& calB,
      const strip_spacepoint_calibration& calT, const lin_circle& lt,
      const seedfinder_config& config, scalar& curvature,
      scalar& impact_parameter);

  /// Check if two doublets with common middle spacepoint can form a triplet
  ///
  /// @param spM is middle spacepoint
  /// @param lb is transformed coordinate of middle-bottom doublet
  /// @param lt is transformed coordinate of middle-top doublet
  /// @param config is configuration parameter
  /// @param iSinTheta2 is the square of sin of pitch angle
  /// @param scatteringInRegion2 is the threshold for scattering angle for the
  /// lower pT cut
  /// @param curvature is curvature of triplet
  /// @param impact_parameter is impact parameter of triplet
  ///
  /// @return boolean value for compatibility
  template <typename T>
  static inline TRACCC_HOST_DEVICE bool isCompatible(
      const edm::spacepoint<T>& spM, const lin_circle& lb, const lin_circle& lt,
      const seedfinder_config& config, const scalar& iSinTheta2,
      const scalar& scatteringInRegion2, scalar& curvature,
      scalar& impact_parameter);
};

template <typename T>
bool TRACCC_HOST_DEVICE triplet_finding_helper::isCompatible(
    const edm::spacepoint<T>& spM, const lin_circle& lb, const lin_circle& lt,
    const seedfinder_config& config, const scalar& iSinTheta2,
    const scalar& scatteringInRegion2, scalar& curvature,
    scalar& impact_parameter) {
  // add errors of spB-spM and spM-spT pairs and add the correlation term
  // for errors on spM
  scalar error2 = lt.Er() + lb.Er() +
                  static_cast<scalar>(2.f) *
                      (lb.cotTheta() * lt.cotTheta() * spM.radius_variance() +
                       spM.z_variance()) *
                      lb.iDeltaR() * lt.iDeltaR();

  scalar deltaCotTheta = lb.cotTheta() - lt.cotTheta();
  scalar deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
  scalar error{0.f};
  scalar dCotThetaMinusError2{0.f};

  // if the error is larger than the difference in theta, no need to
  // compare with scattering
  if (deltaCotTheta2 - error2 > 0) {
    deltaCotTheta = math::fabs(deltaCotTheta);
    // if deltaTheta larger than the scattering for the lower pT cut, skip
    error = std::sqrt(error2);
    dCotThetaMinusError2 = deltaCotTheta2 + error2 -
                           static_cast<scalar>(2.) * deltaCotTheta * error;
    // avoid taking root of scatteringInRegion
    // if left side of ">" is positive, both sides of inequality can be
    // squared
    // (scattering is always positive)
    if (dCotThetaMinusError2 > scatteringInRegion2) {
      return false;
    }
  }

  // protects against division by 0
  scalar dU = lt.U() - lb.U();
  if (dU == static_cast<scalar>(0.)) {
    return false;
  }

  // A and B are evaluated as a function of the circumference parameters
  // x_0 and y_0
  scalar A = (lt.V() - lb.V()) / dU;
  scalar S2 = static_cast<scalar>(1.) + A * A;
  scalar B = lb.V() - A * lb.U();
  scalar B2 = B * B;
  // sqrt(S2)/B = 2 * helixradius
  // calculated radius must not be smaller than minimum radius
  if (S2 < B2 * config.minHelixDiameter2) {
    return false;
  }

  // 1/helixradius: (B/sqrt(S2))*2 (we leave everything squared)
  scalar iHelixDiameter2 = B2 / S2;
  // calculate scattering for p(T) calculated from seed curvature
  scalar pT2scatter =
      static_cast<scalar>(4.) * iHelixDiameter2 * config.pT2perRadius;
  // if pT > maxPtScattering, calculate allowed scattering angle using
  // maxPtScattering instead of pt.
  scalar pT =
      config.pTPerHelixRadius * std::sqrt(S2 / B2) / static_cast<scalar>(2.);
  if (pT > config.maxPtScattering) {
    scalar pTscatter = config.highland / config.maxPtScattering;
    pT2scatter = pTscatter * pTscatter;
  }
  // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
  // from rad to deltaCotTheta
  scalar p2scatter = pT2scatter * iSinTheta2;
  // if deltaTheta larger than allowed scattering for calculated pT, skip
  if ((deltaCotTheta2 - error2 > static_cast<scalar>(0.)) &&
      (dCotThetaMinusError2 >
       p2scatter * config.sigmaScattering * config.sigmaScattering)) {
    return false;
  }

  // calculate curvature
  curvature = B / std::sqrt(S2);

  // A and B allow calculation of impact params in U/V plane with linear
  // function
  // (in contrast to having to solve a quadratic function in x/y plane)
  impact_parameter = math::fabs((A - B * spM.radius()) * spM.radius());

  if (impact_parameter > config.impactMax) {
    return false;
  }

  return true;
}

inline bool TRACCC_HOST_DEVICE triplet_finding_helper::isCompatibleStrip(
    const strip_triplet_context& ctx, const strip_spacepoint_calibration& calM,
    const strip_spacepoint_calibration& calB,
    const strip_spacepoint_calibration& calT, const lin_circle& lt,
    const seedfinder_config& config, scalar& curvature,
    scalar& impact_parameter) {
  // protects against division by 0
  const scalar dU0 = lt.U() - ctx.Ub0;
  if (dU0 == 0.f) {
    return false;
  }
  // A and B are evaluated as a function of the circumference parameters
  // x_0 and y_0
  const scalar A0 = (lt.V() - ctx.Vb0) / dU0;

  // The middle strip check is scale-invariant, so cosTheta is used as the
  // z-component instead of cosTheta * sqrt(1 + A0^2)
  const vector3 directionMiddle = {
      ctx.rotationTermsUVtoXY[0] - ctx.rotationTermsUVtoXY[1] * A0,
      ctx.rotationTermsUVtoXY[0] * A0 + ctx.rotationTermsUVtoXY[1],
      ctx.cosTheta};
  point3 rMTransf{};
  if (!calibrate_strip_spacepoint(directionMiddle, calM, rMTransf,
                                  config.toleranceParam)) {
    return false;
  }
  const scalar zDirectionMiddle = ctx.cosTheta * std::sqrt(1.f + A0 * A0);

  // coordinate transformation and checks for bottom space point
  const scalar B0 = 2.f * (ctx.Vb0 - A0 * ctx.Ub0);
  const scalar Cb = 1.f - B0 * ctx.yB0;
  const scalar Sb = A0 + B0 * ctx.xB0;
  const vector3 directionBottom = {
      ctx.rotationTermsUVtoXY[0] * Cb - ctx.rotationTermsUVtoXY[1] * Sb,
      ctx.rotationTermsUVtoXY[0] * Sb + ctx.rotationTermsUVtoXY[1] * Cb,
      zDirectionMiddle};
  point3 rBTransf{};
  if (!calibrate_strip_spacepoint(directionBottom, calB, rBTransf,
                                  config.toleranceParam)) {
    return false;
  }

  // coordinate transformation and checks for top space point
  const scalar xT0 = lt.U() / (lt.iDeltaR() * lt.iDeltaR());
  const scalar yT0 = lt.V() / (lt.iDeltaR() * lt.iDeltaR());
  const scalar Ct = 1.f - B0 * yT0;
  const scalar St = A0 + B0 * xT0;
  const vector3 directionTop = {
      ctx.rotationTermsUVtoXY[0] * Ct - ctx.rotationTermsUVtoXY[1] * St,
      ctx.rotationTermsUVtoXY[0] * St + ctx.rotationTermsUVtoXY[1] * Ct,
      zDirectionMiddle};
  point3 rTTransf{};
  if (!calibrate_strip_spacepoint(directionTop, calT, rTTransf,
                                  config.toleranceParam)) {
    return false;
  }

  // bottom and top coordinates in the spM reference frame
  const scalar xB = rBTransf[0] - rMTransf[0];
  const scalar yB = rBTransf[1] - rMTransf[1];
  const scalar zB = rBTransf[2] - rMTransf[2];
  const scalar xT = rTTransf[0] - rMTransf[0];
  const scalar yT = rTTransf[1] - rMTransf[1];
  const scalar zT = rTTransf[2] - rMTransf[2];

  const scalar iDeltaRB2 = 1.f / (xB * xB + yB * yB);
  const scalar iDeltaRT2 = 1.f / (xT * xT + yT * yT);

  const scalar cotThetaB = -zB * std::sqrt(iDeltaRB2);
  const scalar cotThetaT = zT * std::sqrt(iDeltaRT2);

  // use arithmetic average
  const scalar averageCotTheta = 0.5f * (cotThetaB + cotThetaT);
  const scalar cotThetaAvg2 = averageCotTheta * averageCotTheta;

  // add errors of spB-spM and spM-spT pairs and add the correlation term
  // for errors on spM
  const scalar error2 = lt.Er() + ctx.erB +
                        2.f * (cotThetaAvg2 * ctx.varianceRM + ctx.varianceZM) *
                            ctx.iDeltaRB * lt.iDeltaR();

  const scalar deltaCotTheta = cotThetaB - cotThetaT;
  const scalar deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

  // compatibility of the r-z slopes of the two seed segments within the
  // measurement errors and the scattering for the lower pT cut
  if (deltaCotTheta2 > error2 + ctx.scatteringInRegion2) {
    return false;
  }

  const scalar rMxy =
      std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
  const scalar irMxy = 1.f / rMxy;
  const scalar Ax = rMTransf[0] * irMxy;
  const scalar Ay = rMTransf[1] * irMxy;

  const scalar Ub = (xB * Ax + yB * Ay) * iDeltaRB2;
  const scalar Vb = (yB * Ax - xB * Ay) * iDeltaRB2;
  const scalar Ut = (xT * Ax + yT * Ay) * iDeltaRT2;
  const scalar Vt = (yT * Ax - xT * Ay) * iDeltaRT2;

  const scalar dU = Ut - Ub;
  // protects against division by 0
  if (dU == 0.f) {
    return false;
  }
  const scalar A = (Vt - Vb) / dU;
  const scalar S2 = 1.f + A * A;
  const scalar B = Vb - A * Ub;
  const scalar B2 = B * B;

  // sqrt(S2)/B = 2 * helixradius
  // calculated radius must not be smaller than minimum radius
  if (S2 < B2 * config.minHelixDiameter2) {
    return false;
  }

  // refinement of the slope compatibility using the scattering for the
  // measured pT
  const scalar iHelixDiameter2 = B2 / S2;
  const scalar p2scatterSigma = iHelixDiameter2 * ctx.sigmaSquaredPtDependent;
  if (deltaCotTheta2 > error2 + p2scatterSigma) {
    return false;
  }

  curvature = B / std::sqrt(S2);

  // A and B allow calculation of impact params in U/V plane with linear
  // function
  impact_parameter = math::fabs((A - B * rMxy) * rMxy);
  if (impact_parameter > config.impactMax) {
    return false;
  }

  return true;
}

/// Make the triplet compatibility check for the top candidates of one
/// middle-bottom doublet
///
/// The check is chosen at compile time, so that the strip aware one, and the
/// doublet quantities it needs, only exist in code that uses them.
///
/// @tparam use_strip_info whether to apply the strip aware check
/// @param spM is the middle spacepoint
/// @param lb is transformed coordinate of middle-bottom doublet
/// @param strip_calibrations are the strip calibrations of all spacepoints,
///        only read by the strip aware check
/// @param spM_idx is the index of the middle spacepoint
/// @param spB_idx is the index of the bottom spacepoint
/// @param config is configuration parameter
/// @param iSinTheta2 is the inverse of sin^2(theta) of the middle-bottom
///        doublet
/// @param scatteringInRegion2 is the scattering term for the lower pT cut
///
/// @return a callable taking the top spacepoint index, the transformed
///         coordinate of the middle-top doublet, and the curvature and impact
///         parameter outputs
template <bool use_strip_info, typename T>
TRACCC_HOST_DEVICE inline auto make_triplet_compatibility_check(
    const edm::spacepoint<T>& spM, const lin_circle& lb,
    const strip_spacepoint_calibration_collection_types::const_device&
        strip_calibrations,
    const unsigned int spM_idx, const unsigned int spB_idx,
    const seedfinder_config& config, const scalar iSinTheta2,
    const scalar scatteringInRegion2) {
  if constexpr (use_strip_info) {
    return [ctx = strip_triplet_context(spM, lb, config, iSinTheta2,
                                        scatteringInRegion2),
            &calM = strip_calibrations.at(spM_idx),
            &calB = strip_calibrations.at(spB_idx), &strip_calibrations,
            &config](const unsigned int spT_idx, const lin_circle& lt,
                     scalar& curvature, scalar& impact_parameter) {
      return triplet_finding_helper::isCompatibleStrip(
          ctx, calM, calB, strip_calibrations.at(spT_idx), lt, config,
          curvature, impact_parameter);
    };
  } else {
    return [&spM, &lb, &config, iSinTheta2, scatteringInRegion2](
               const unsigned int, const lin_circle& lt, scalar& curvature,
               scalar& impact_parameter) {
      return triplet_finding_helper::isCompatible(
          spM, lb, lt, config, iSinTheta2, scatteringInRegion2, curvature,
          impact_parameter);
    };
  }
}

}  // namespace traccc
