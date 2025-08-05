// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <Eigen/Dense>
#include <boost/mp11.hpp>
#include <boost/mp11/algorithm.hpp>

namespace Acts::Experimental {

namespace {

/// Check the compatibility of strip space point coordinates in xyz assuming
/// the Bottom-Middle direction with the strip measurement details
static bool stripCoordinateCheck(float tolerance,
                                 const ConstSpacePointProxy2& sp,
                                 const Eigen::Vector3f& spacePointPosition,
                                 Eigen::Vector3f& outputCoordinates) {
  const Eigen::Vector3f& topStripVector = sp.topStripVector();
  const Eigen::Vector3f& bottomStripVector = sp.bottomStripVector();
  const Eigen::Vector3f& stripCenterDistance = sp.stripCenterDistance();

  // cross product between top strip vector and spacepointPosition
  Eigen::Vector3f d1 = topStripVector.cross(spacePointPosition);

  // scalar product between bottom strip vector and d1
  float bd1 = bottomStripVector.dot(d1);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  float s1 = stripCenterDistance.dot(d1);
  if (std::abs(s1) > std::abs(bd1) * tolerance) {
    return false;
  }

  // cross product between bottom strip vector and spacepointPosition
  Eigen::Vector3f d0 = bottomStripVector.cross(spacePointPosition);

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the top detector element
  float s0 = stripCenterDistance.dot(d0);
  if (std::abs(s0) > std::abs(bd1) * tolerance) {
    return false;
  }

  // if arrive here spacepointPosition is compatible with strip directions and
  // detector elements

  const Eigen::Vector3f& topStripCenter = sp.topStripCenter();

  // spacepointPosition corrected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates = topStripCenter + topStripVector * s0;
  return true;
}

}  // namespace

template <bool useStripInfo, bool sortedInCotTheta>
class TripletSeedFinder::Impl final : public TripletSeedFinder::ImplBase {
 public:
  using ImplBase::ImplBase;

  template <typename TopDoublets>
  void createPixelTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet, TopDoublets& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    const float rM = spM.r();
    const float varianceRM = spM.varianceR();
    const float varianceZM = spM.varianceZ();

    const LinCircle& lb = bottomDoublet.linCircle();

    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1 + cotThetaB * cotThetaB;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    float scatteringInRegion2 = m_cfg.multipleScattering2 * iSinTheta2;

    // Reserve enough space, in case current capacity is too little
    tripletTopCandidates.reserve(tripletTopCandidates.size() +
                                 topDoublets.size());

    for (std::size_t topDoubletIndex = 0; topDoubletIndex < topDoublets.size();
         ++topDoubletIndex) {
      auto topDoublet = topDoublets[topDoubletIndex];
      const ConstSpacePointProxy2 spT = spacePoints[topDoublet.spacePoint()];
      const LinCircle& lt = topDoublet.linCircle();
      float cotThetaT = lt.cotTheta;

      // use geometric average
      float cotThetaAvg2 = cotThetaB * cotThetaT;
      if (cotThetaAvg2 <= 0) {
        continue;
      }

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 =
          lt.Er + ErB +
          2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB * lt.iDeltaR;

      float deltaCotTheta = cotThetaB - cotThetaT;
      float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

      // Apply a cut on the compatibility between the r-z slope of the two
      // seed segments. This is done by comparing the squared difference
      // between slopes, and comparing to the squared uncertainty in this
      // difference - we keep a seed if the difference is compatible within
      // the assumed uncertainties. The uncertainties get contribution from
      // the  space-point-related squared error (error2) and a scattering term
      // calculated assuming the minimum pt we expect to reconstruct
      // (scatteringInRegion2). This assumes gaussian error propagation which
      // allows just adding the two errors if they are uncorrelated (which is
      // fair for scattering and measurement uncertainties)
      if (deltaCotTheta2 > error2 + scatteringInRegion2) {
        if constexpr (sortedInCotTheta) {
          // skip top SPs based on cotTheta sorting when producing triplets
          // break if cotTheta from bottom SP < cotTheta from top SP because
          // the SP are sorted by cotTheta
          if (cotThetaB < cotThetaT) {
            break;
          }
          topDoublets = topDoublets.subrange(topDoubletIndex + 1);
        }
        continue;
      }

      float dU = lt.U - Ub;
      // protects against division by 0
      if (dU == 0.) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A = (lt.V - Vb) / dU;
      float S2 = 1 + A * A;
      float B = Vb - A * Ub;
      float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_cfg.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      float sigmaSquaredPtDependent = iSinTheta2 * m_cfg.sigmapT2perRadius;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      if (!std::isinf(m_cfg.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 || m_cfg.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * m_cfg.maxPtScattering) {
          float pTscatterSigma =
              (m_cfg.highland / m_cfg.maxPtScattering) * m_cfg.sigmaScattering;
          p2scatterSigma = pTscatterSigma * pTscatterSigma * iSinTheta2;
        }
      }

      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > error2 + p2scatterSigma) {
        if constexpr (sortedInCotTheta) {
          if (cotThetaB < cotThetaT) {
            break;
          }
          topDoublets = topDoublets.subrange(topDoubletIndex);
        }
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rM) * rM);
      if (Im > m_cfg.impactMax) {
        continue;
      }

      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      tripletTopCandidates.emplace_back(spT.index(), B / std::sqrt(S2), Im);
    }  // loop on tops
  }

  template <typename TopDoublets>
  void createStripTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet, TopDoublets& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    float rM = spM.r();
    float cosPhiM = spM.x() / rM;
    float sinPhiM = spM.y() / rM;
    float varianceRM = spM.varianceR();
    float varianceZM = spM.varianceZ();

    // Reserve enough space, in case current capacity is too little
    tripletTopCandidates.reserve(tripletTopCandidates.size() +
                                 topDoublets.size());

    const ConstSpacePointProxy2 spB = spacePoints[bottomDoublet.spacePoint()];
    const LinCircle& lb = bottomDoublet.linCircle();

    float cotThetaB = lb.cotTheta;
    float Vb = lb.V;
    float Ub = lb.U;
    float ErB = lb.Er;
    float iDeltaRB = lb.iDeltaR;

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    float iSinTheta2 = 1 + cotThetaB * cotThetaB;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    float scatteringInRegion2 = m_cfg.multipleScattering2 * iSinTheta2;

    float sinTheta = 1 / std::sqrt(iSinTheta2);
    float cosTheta = cotThetaB * sinTheta;

    // coordinate transformation and checks for middle spacepoint
    // x and y terms for the rotation from UV to XY plane
    Eigen::Vector2f rotationTermsUVtoXY = {cosPhiM * sinTheta,
                                           sinPhiM * sinTheta};

    for (auto topDoublet : topDoublets) {
      const ConstSpacePointProxy2 spT = spacePoints[topDoublet.spacePoint()];
      const LinCircle& lt = topDoublet.linCircle();

      // protects against division by 0
      float dU = lt.U - Ub;
      if (dU == 0.) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A0 = (lt.V - Vb) / dU;

      float zPositionMiddle = cosTheta * std::sqrt(1 + A0 * A0);

      // position of Middle SP converted from UV to XY assuming cotTheta
      // evaluated from the Bottom and Middle SPs double
      Eigen::Vector3f positionMiddle = {
          rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
          rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
          zPositionMiddle};

      Eigen::Vector3f rMTransf;
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spM, positionMiddle,
                                rMTransf)) {
        continue;
      }

      // coordinate transformation and checks for bottom spacepoint
      float B0 = 2 * (Vb - A0 * Ub);
      float Cb = 1 - B0 * lb.y;
      float Sb = A0 + B0 * lb.x;
      Eigen::Vector3f positionBottom = {
          rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
          rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
          zPositionMiddle};

      Eigen::Vector3f rBTransf;
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spB, positionBottom,
                                rBTransf)) {
        continue;
      }

      // coordinate transformation and checks for top spacepoint
      float Ct = 1 - B0 * lt.y;
      float St = A0 + B0 * lt.x;
      Eigen::Vector3f positionTop = {
          rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
          rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
          zPositionMiddle};

      Eigen::Vector3f rTTransf;
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spT, positionTop,
                                rTTransf)) {
        continue;
      }

      // bottom and top coordinates in the spM reference frame
      float xB = rBTransf[0] - rMTransf[0];
      float yB = rBTransf[1] - rMTransf[1];
      float zB = rBTransf[2] - rMTransf[2];
      float xT = rTTransf[0] - rMTransf[0];
      float yT = rTTransf[1] - rMTransf[1];
      float zT = rTTransf[2] - rMTransf[2];

      float iDeltaRB2 = 1 / (xB * xB + yB * yB);
      float iDeltaRT2 = 1 / (xT * xT + yT * yT);

      cotThetaB = -zB * std::sqrt(iDeltaRB2);
      float cotThetaT = zT * std::sqrt(iDeltaRT2);

      // use arithmetic average
      float averageCotTheta = 0.5f * (cotThetaB + cotThetaT);
      float cotThetaAvg2 = averageCotTheta * averageCotTheta;

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 =
          lt.Er + ErB +
          2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB * lt.iDeltaR;

      float deltaCotTheta = cotThetaB - cotThetaT;
      float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

      // Apply a cut on the compatibility between the r-z slope of the two
      // seed segments. This is done by comparing the squared difference
      // between slopes, and comparing to the squared uncertainty in this
      // difference - we keep a seed if the difference is compatible within
      // the assumed uncertainties. The uncertainties get contribution from
      // the  space-point-related squared error (error2) and a scattering term
      // calculated assuming the minimum pt we expect to reconstruct
      // (scatteringInRegion2). This assumes gaussian error propagation which
      // allows just adding the two errors if they are uncorrelated (which is
      // fair for scattering and measurement uncertainties)
      if (deltaCotTheta2 > error2 + scatteringInRegion2) {
        // skip top SPs based on cotTheta sorting when producing triplets
        continue;
      }

      float rMxy =
          std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
      float irMxy = 1 / rMxy;
      float Ax = rMTransf[0] * irMxy;
      float Ay = rMTransf[1] * irMxy;

      float ub = (xB * Ax + yB * Ay) * iDeltaRB2;
      float vb = (yB * Ax - xB * Ay) * iDeltaRB2;
      float ut = (xT * Ax + yT * Ay) * iDeltaRT2;
      float vt = (yT * Ax - xT * Ay) * iDeltaRT2;

      dU = ut - ub;
      // protects against division by 0
      if (dU == 0.) {
        continue;
      }
      float A = (vt - vb) / dU;
      float S2 = 1 + A * A;
      float B = vb - A * ub;
      float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_cfg.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      float iHelixDiameter2 = B2 / S2;
      float sigmaSquaredPtDependent = iSinTheta2 * m_cfg.sigmapT2perRadius;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      if (!std::isinf(m_cfg.maxPtScattering)) {
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        // To avoid 0-divison the pT check is skipped in case of B2==0, and
        // p2scatterSigma is calculated directly from maxPtScattering
        if (B2 == 0 || m_cfg.pTPerHelixRadius * std::sqrt(S2 / B2) >
                           2. * m_cfg.maxPtScattering) {
          float pTscatterSigma =
              (m_cfg.highland / m_cfg.maxPtScattering) * m_cfg.sigmaScattering;
          p2scatterSigma = pTscatterSigma * pTscatterSigma * iSinTheta2;
        }
      }

      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > error2 + p2scatterSigma) {
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rMxy) * rMxy);
      if (Im > m_cfg.impactMax) {
        continue;
      }

      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      tripletTopCandidates.emplace_back(spT.index(), B / std::sqrt(S2), Im);
    }  // loop on tops
  }

  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Range& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const override {
    if constexpr (useStripInfo) {
      createStripTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    } else {
      createPixelTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    }
  }

  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const override {
    if constexpr (useStripInfo) {
      createStripTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    } else {
      createPixelTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    }
  }

  void createTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      DoubletsForMiddleSp::Subset2& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const override {
    if constexpr (useStripInfo) {
      createStripTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    } else {
      createPixelTripletTopCandidates(spacePoints, spM, bottomDoublet,
                                      topDoublets, tripletTopCandidates);
    }
  }
};

std::shared_ptr<TripletSeedFinder::ImplBase> TripletSeedFinder::makeImpl(
    const DerivedConfig& config) {
  using BooleanOptions =
      boost::mp11::mp_list<std::bool_constant<false>, std::bool_constant<true>>;

  using UseStripInfoOptions = BooleanOptions;
  using SortedInCotThetaOptions = BooleanOptions;

  using TripletOptions =
      boost::mp11::mp_product<boost::mp11::mp_list, UseStripInfoOptions,
                              SortedInCotThetaOptions>;

  std::shared_ptr<TripletSeedFinder::ImplBase> result;
  boost::mp11::mp_for_each<TripletOptions>([&](auto option) {
    using OptionType = decltype(option);

    using UseStripInfo = boost::mp11::mp_at_c<OptionType, 0>;
    using SortedInCotTheta = boost::mp11::mp_at_c<OptionType, 1>;

    if (config.useStripInfo != UseStripInfo::value ||
        config.sortedByCotTheta != SortedInCotTheta::value) {
      return;  // skip if the configuration does not match
    }

    // check if we already have an implementation for this configuration
    if (result != nullptr) {
      throw std::runtime_error(
          "TripletSeedFinder: Multiple implementations found for one "
          "configuration");
    }

    // create the implementation for the given configuration
    result = std::make_shared<
        TripletSeedFinder::Impl<UseStripInfo::value, SortedInCotTheta::value>>(
        config);
  });
  if (result == nullptr) {
    throw std::runtime_error(
        "TripletSeedFinder: No implementation found for the given "
        "configuration");
  }
  return result;
}

TripletSeedFinder::DerivedConfig::DerivedConfig(const Config& config,
                                                float bFieldInZ_)
    : Config(config), bFieldInZ(bFieldInZ_) {
  using namespace Acts::UnitLiterals;

  // similar to `theta0Highland` in `Core/src/Material/Interactions.cpp`
  {
    const double xOverX0 = radLengthPerSeed;
    const double q2OverBeta2 = 1;  // q^2=1, beta^2~1
    // RPP2018 eq. 33.15 (treats beta and q² consistently)
    const double t = std::sqrt(xOverX0 * q2OverBeta2);
    // log((x/X0) * (q²/beta²)) = log((sqrt(x/X0) * (q/beta))²)
    //                          = 2 * log(sqrt(x/X0) * (q/beta))
    highland =
        static_cast<float>(13.6_MeV * t * (1.0 + 0.038 * 2 * std::log(t)));
  }

  const float maxScatteringAngle = highland / minPt;
  const float maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // bFieldInZ is in (pT/radius) natively, no need for conversion
  pTPerHelixRadius = bFieldInZ;
  minHelixDiameter2 = square(minPt * 2 / pTPerHelixRadius) * helixCutTolerance;
  const float pT2perRadius = square(highland / pTPerHelixRadius);
  sigmapT2perRadius = pT2perRadius * square(2 * sigmaScattering);
  multipleScattering2 = maxScatteringAngle2 * square(sigmaScattering);
}

TripletSeedFinder::TripletSeedFinder(const DerivedConfig& config)
    : m_impl(makeImpl(config)) {}

}  // namespace Acts::Experimental
