// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/TripletSeedFinder.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <ranges>

#include <Eigen/Dense>
#include <boost/mp11.hpp>
#include <boost/mp11/algorithm.hpp>

namespace Acts {

namespace {

/// Check the compatibility of strip space point coordinates in xyz assuming
/// the Bottom-Middle direction with the strip measurement details
bool stripCoordinateCheck(float tolerance, const ConstSpacePointProxy2& sp,
                          const std::array<float, 3>& spacePointPosition,
                          std::array<float, 3>& outputCoordinates) {
  const std::array<float, 3>& topStripVector = sp.topStripVector();
  const std::array<float, 3>& bottomStripVector = sp.bottomStripVector();
  const std::array<float, 3>& stripCenterDistance = sp.stripCenterDistance();

  // cross product between top strip vector and spacePointPosition
  const std::array<float, 3> d1 = {
      topStripVector[1] * spacePointPosition[2] -
          topStripVector[2] * spacePointPosition[1],
      topStripVector[2] * spacePointPosition[0] -
          topStripVector[0] * spacePointPosition[2],
      topStripVector[0] * spacePointPosition[1] -
          topStripVector[1] * spacePointPosition[0]};

  // scalar product between bottom strip vector and d1
  const float bd1 = bottomStripVector[0] * d1[0] +
                    bottomStripVector[1] * d1[1] + bottomStripVector[2] * d1[2];

  // compatibility check using distance between strips to evaluate if
  // spacepointPosition is inside the bottom detector element
  const float s1 = stripCenterDistance[0] * d1[0] +
                   stripCenterDistance[1] * d1[1] +
                   stripCenterDistance[2] * d1[2];
  if (std::abs(s1) > std::abs(bd1) * tolerance) {
    return false;
  }

  // cross product between bottom strip vector and spacePointPosition
  const std::array<float, 3> d0 = {
      bottomStripVector[1] * spacePointPosition[2] -
          bottomStripVector[2] * spacePointPosition[1],
      bottomStripVector[2] * spacePointPosition[0] -
          bottomStripVector[0] * spacePointPosition[2],
      bottomStripVector[0] * spacePointPosition[1] -
          bottomStripVector[1] * spacePointPosition[0]};

  // compatibility check using distance between strips to evaluate if
  // spacePointPosition is inside the top detector element
  float s0 = stripCenterDistance[0] * d0[0] + stripCenterDistance[1] * d0[1] +
             stripCenterDistance[2] * d0[2];
  if (std::abs(s0) > std::abs(bd1) * tolerance) {
    return false;
  }

  // if arrive here spacePointPosition is compatible with strip directions and
  // detector elements

  const std::array<float, 3>& topStripCenter = sp.topStripCenter();

  // spacePointPosition corrected with respect to the top strip position and
  // direction and the distance between the strips
  s0 = s0 / bd1;
  outputCoordinates[0] = topStripCenter[0] + topStripVector[0] * s0;
  outputCoordinates[1] = topStripCenter[1] + topStripVector[1] * s0;
  outputCoordinates[2] = topStripCenter[2] + topStripVector[2] * s0;
  return true;
}

template <bool useStripInfo, bool sortedByCotTheta>
class Impl final : public TripletSeedFinder {
 public:
  explicit Impl(const DerivedConfig& config) : m_cfg(config) {}

  const DerivedConfig& config() const override { return m_cfg; }

  template <typename TopDoublets>
  void createPixelTripletTopCandidates(
      const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet, TopDoublets& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    const float rM = spM.zr()[1];
    const float varianceZM = spM.varianceZ();
    const float varianceRM = spM.varianceR();

    // Reserve enough space, in case current capacity is too little
    tripletTopCandidates.reserve(tripletTopCandidates.size() +
                                 topDoublets.size());

    const float cotThetaB = bottomDoublet.cotTheta();
    const float erB = bottomDoublet.er();
    const float iDeltaRB = bottomDoublet.iDeltaR();
    const float Ub = bottomDoublet.u();
    const float Vb = bottomDoublet.v();

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    const float iSinTheta2 = 1 + cotThetaB * cotThetaB;
    const float sigmaSquaredPtDependent = iSinTheta2 * m_cfg.sigmapT2perRadius;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    const float scatteringInRegion2 = m_cfg.multipleScattering2 * iSinTheta2;

    std::size_t topDoubletOffset = 0;
    for (auto [topDoublet, topDoubletIndex] :
         zip(topDoublets, std::ranges::iota_view<std::size_t, std::size_t>(
                              0, topDoublets.size()))) {
      const SpacePointIndex2 spT = topDoublet.spacePointIndex();
      const float cotThetaT = topDoublet.cotTheta();

      // use geometric average
      const float cotThetaAvg2 = cotThetaB * cotThetaT;

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      const float error2 = topDoublet.er() + erB +
                           2 * (cotThetaAvg2 * varianceRM + varianceZM) *
                               iDeltaRB * topDoublet.iDeltaR();

      const float deltaCotTheta = cotThetaB - cotThetaT;
      const float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

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
        if constexpr (sortedByCotTheta) {
          // skip top SPs based on cotTheta sorting when producing triplets
          // break if cotTheta from bottom SP < cotTheta from top SP because
          // the SP are sorted by cotTheta
          if (cotThetaB < cotThetaT) {
            break;
          }
          topDoubletOffset = topDoubletIndex + 1;
        }
        continue;
      }

      const float dU = topDoublet.u() - Ub;
      // protects against division by 0
      if (dU == 0) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      const float A = (topDoublet.v() - Vb) / dU;
      const float S2 = 1 + A * A;
      const float B = Vb - A * Ub;
      const float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_cfg.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      const float iHelixDiameter2 = B2 / S2;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      const float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > error2 + p2scatterSigma) {
        if constexpr (sortedByCotTheta) {
          if (cotThetaB < cotThetaT) {
            break;
          }
          topDoubletOffset = topDoubletIndex;
        }
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      const float im = std::abs((A - B * rM) * rM);
      if (im > m_cfg.impactMax) {
        continue;
      }

      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      tripletTopCandidates.emplace_back(spT, B / std::sqrt(S2), im);
    }  // loop on tops

    if constexpr (sortedByCotTheta) {
      // remove the top doublets that were skipped due to cotTheta sorting
      topDoublets = topDoublets.subrange(topDoubletOffset);
    }
  }

  template <typename TopDoublets>
  void createStripTripletTopCandidates(
      const SpacePointContainer2& spacePoints, const ConstSpacePointProxy2& spM,
      const DoubletsForMiddleSp::Proxy& bottomDoublet,
      const TopDoublets& topDoublets,
      TripletTopCandidates& tripletTopCandidates) const {
    const float rM = spM.zr()[1];
    const float cosPhiM = spM.xy()[0] / rM;
    const float sinPhiM = spM.xy()[1] / rM;
    const float varianceZM = spM.varianceZ();
    const float varianceRM = spM.varianceR();

    // Reserve enough space, in case current capacity is too little
    tripletTopCandidates.reserve(tripletTopCandidates.size() +
                                 topDoublets.size());

    float cotThetaB = bottomDoublet.cotTheta();
    const float erB = bottomDoublet.er();
    const float iDeltaRB = bottomDoublet.iDeltaR();
    const float Vb = bottomDoublet.v();
    const float Ub = bottomDoublet.u();

    // 1+(cot^2(theta)) = 1/sin^2(theta)
    const float iSinTheta2 = 1 + cotThetaB * cotThetaB;
    const float sigmaSquaredPtDependent = iSinTheta2 * m_cfg.sigmapT2perRadius;
    // calculate max scattering for min momentum at the seed's theta angle
    // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
    // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
    // scattering
    // but to avoid trig functions we approximate cot by scaling by
    // 1/sin^4(theta)
    // resolving with pT to p scaling --> only divide by sin^2(theta)
    // max approximation error for allowed scattering angles of 0.04 rad at
    // eta=infinity: ~8.5%
    const float scatteringInRegion2 = m_cfg.multipleScattering2 * iSinTheta2;

    const float sinTheta = 1 / std::sqrt(iSinTheta2);
    const float cosTheta = cotThetaB * sinTheta;

    // coordinate transformation and checks for middle spacepoint
    // x and y terms for the rotation from UV to XY plane
    const std::array<float, 2> rotationTermsUVtoXY = {cosPhiM * sinTheta,
                                                      sinPhiM * sinTheta};

    for (auto topDoublet : topDoublets) {
      // protects against division by 0
      float dU = topDoublet.u() - Ub;
      if (dU == 0) {
        continue;
      }
      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      const float A0 = (topDoublet.v() - Vb) / dU;

      const float zPositionMiddle = cosTheta * std::sqrt(1 + A0 * A0);

      // position of Middle SP converted from UV to XY assuming cotTheta
      // evaluated from the Bottom and Middle SPs double
      const std::array<float, 3> positionMiddle = {
          rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
          rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
          zPositionMiddle};

      std::array<float, 3> rMTransf{};
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spM, positionMiddle,
                                rMTransf)) {
        continue;
      }

      // coordinate transformation and checks for bottom spacepoint
      const float B0 = 2 * (Vb - A0 * Ub);
      const float Cb = 1 - B0 * bottomDoublet.y();
      const float Sb = A0 + B0 * bottomDoublet.x();
      const std::array<float, 3> positionBottom = {
          rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
          rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
          zPositionMiddle};

      const ConstSpacePointProxy2 spB =
          spacePoints[bottomDoublet.spacePointIndex()];
      std::array<float, 3> rBTransf{};
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spB, positionBottom,
                                rBTransf)) {
        continue;
      }

      // coordinate transformation and checks for top spacepoint
      const float Ct = 1 - B0 * topDoublet.y();
      const float St = A0 + B0 * topDoublet.x();
      const std::array<float, 3> positionTop = {
          rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
          rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
          zPositionMiddle};

      const ConstSpacePointProxy2 spT =
          spacePoints[topDoublet.spacePointIndex()];
      std::array<float, 3> rTTransf{};
      if (!stripCoordinateCheck(m_cfg.toleranceParam, spT, positionTop,
                                rTTransf)) {
        continue;
      }

      // bottom and top coordinates in the spM reference frame
      const float xB = rBTransf[0] - rMTransf[0];
      const float yB = rBTransf[1] - rMTransf[1];
      const float zB = rBTransf[2] - rMTransf[2];
      const float xT = rTTransf[0] - rMTransf[0];
      const float yT = rTTransf[1] - rMTransf[1];
      const float zT = rTTransf[2] - rMTransf[2];

      const float iDeltaRB2 = 1 / (xB * xB + yB * yB);
      const float iDeltaRT2 = 1 / (xT * xT + yT * yT);

      cotThetaB = -zB * std::sqrt(iDeltaRB2);
      const float cotThetaT = zT * std::sqrt(iDeltaRT2);

      // use arithmetic average
      const float averageCotTheta = 0.5f * (cotThetaB + cotThetaT);
      const float cotThetaAvg2 = averageCotTheta * averageCotTheta;

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      const float error2 = topDoublet.er() + erB +
                           2 * (cotThetaAvg2 * varianceRM + varianceZM) *
                               iDeltaRB * topDoublet.iDeltaR();

      const float deltaCotTheta = cotThetaB - cotThetaT;
      const float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;

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

      const float rMxy =
          std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
      const float irMxy = 1 / rMxy;
      const float Ax = rMTransf[0] * irMxy;
      const float Ay = rMTransf[1] * irMxy;

      const float ub = (xB * Ax + yB * Ay) * iDeltaRB2;
      const float vb = (yB * Ax - xB * Ay) * iDeltaRB2;
      const float ut = (xT * Ax + yT * Ay) * iDeltaRT2;
      const float vt = (yT * Ax - xT * Ay) * iDeltaRT2;

      dU = ut - ub;
      // protects against division by 0
      if (dU == 0) {
        continue;
      }
      const float A = (vt - vb) / dU;
      const float S2 = 1 + A * A;
      const float B = vb - A * ub;
      const float B2 = B * B;

      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_cfg.minHelixDiameter2) {
        continue;
      }

      // refinement of the cut on the compatibility between the r-z slope of
      // the two seed segments using a scattering term scaled by the actual
      // measured pT (p2scatterSigma)
      const float iHelixDiameter2 = B2 / S2;
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      const float p2scatterSigma = iHelixDiameter2 * sigmaSquaredPtDependent;
      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 > error2 + p2scatterSigma) {
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      const float im = std::abs((A - B * rMxy) * rMxy);
      if (im > m_cfg.impactMax) {
        continue;
      }

      // inverse diameter is signed depending on if the curvature is
      // positive/negative in phi
      tripletTopCandidates.emplace_back(spT.index(), B / std::sqrt(S2), im);
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
      createPixelTripletTopCandidates(spM, bottomDoublet, topDoublets,
                                      tripletTopCandidates);
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
      createPixelTripletTopCandidates(spM, bottomDoublet, topDoublets,
                                      tripletTopCandidates);
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
      createPixelTripletTopCandidates(spM, bottomDoublet, topDoublets,
                                      tripletTopCandidates);
    }
  }

 private:
  DerivedConfig m_cfg;
};

}  // namespace

TripletSeedFinder::DerivedConfig::DerivedConfig(const Config& config,
                                                float bFieldInZ_)
    : Config(config), bFieldInZ(bFieldInZ_) {
  // similar to `theta0Highland` in `Core/src/Material/Interactions.cpp`
  {
    using namespace Acts::UnitLiterals;
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

std::unique_ptr<TripletSeedFinder> TripletSeedFinder::create(
    const DerivedConfig& config) {
  using BooleanOptions =
      boost::mp11::mp_list<std::bool_constant<false>, std::bool_constant<true>>;

  using UseStripInfoOptions = BooleanOptions;
  using SortedByCotThetaOptions = BooleanOptions;

  using TripletOptions =
      boost::mp11::mp_product<boost::mp11::mp_list, UseStripInfoOptions,
                              SortedByCotThetaOptions>;

  std::unique_ptr<TripletSeedFinder> result;
  boost::mp11::mp_for_each<TripletOptions>([&](auto option) {
    using OptionType = decltype(option);

    using UseStripInfo = boost::mp11::mp_at_c<OptionType, 0>;
    using SortedByCotTheta = boost::mp11::mp_at_c<OptionType, 1>;

    if (config.useStripInfo != UseStripInfo::value ||
        config.sortedByCotTheta != SortedByCotTheta::value) {
      return;  // skip if the configuration does not match
    }

    // check if we already have an implementation for this configuration
    if (result != nullptr) {
      throw std::runtime_error(
          "TripletSeedFinder: Multiple implementations found for one "
          "configuration");
    }

    // create the implementation for the given configuration
    result =
        std::make_unique<Impl<UseStripInfo::value, SortedByCotTheta::value>>(
            config);
  });
  if (result == nullptr) {
    throw std::runtime_error(
        "TripletSeedFinder: No implementation found for the given "
        "configuration");
  }
  return result;
}

}  // namespace Acts
