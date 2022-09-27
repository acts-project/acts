// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFilter.hpp"

#include <cmath>
#include <numeric>
#include <type_traits>

#include <chrono>
#include <iostream>

//#include<iostream>

namespace Acts {

template <typename external_spacepoint_t, typename platform_t>
Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(config.toInternalUnits()) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);
  m_config.sigmapT2perRadius =
      m_config.pT2perRadius * std::pow(2 * m_config.sigmaScattering, 2);
}

template <typename external_spacepoint_t, typename platform_t>
template <template <typename...> typename container_t, typename sp_range_t>
void Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    State& state,
    std::back_insert_iterator<container_t<Seed<external_spacepoint_t>>> outIt,
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs,
    const Acts::Vector2 rMiddleSPRange) const {
  for (auto spM : middleSPs) {
    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();

    /// check if spM is outside our radial region of interest
    if (m_config.useVariableMiddleSPRange) {
      if (rM < rMiddleSPRange[0]) {
        continue;
      }
      if (rM > rMiddleSPRange[1]) {
        // break if SP are sorted in r
        if (m_config.forceRadialSorting) {
          break;
        }
        continue;
      }
    } else if (not m_config.rRangeMiddleSP.empty()) {
      /// get zBin position of the middle SP
      auto pVal = std::lower_bound(m_config.zBinEdges.begin(),
                                   m_config.zBinEdges.end(), zM);
      int zBin = std::distance(m_config.zBinEdges.begin(), pVal);
      /// protects against zM at the limit of zBinEdges
      zBin == 0 ? zBin : --zBin;
      if (rM < m_config.rRangeMiddleSP[zBin][0]) {
        continue;
      }
      if (rM > m_config.rRangeMiddleSP[zBin][1]) {
        // break if SP are sorted in r
        if (m_config.forceRadialSorting) {
          break;
        }
        continue;
      }
    }

    size_t nTopSeedConf = 0;
    if (m_config.seedConfirmation == true) {
      // check if middle SP is in the central or forward region
      SeedConfirmationRangeConfig seedConfRange =
          (zM > m_config.centralSeedConfirmationRange.zMaxSeedConf ||
           zM < m_config.centralSeedConfirmationRange.zMinSeedConf)
              ? m_config.forwardSeedConfirmationRange
              : m_config.centralSeedConfirmationRange;
      // set the minimum number of top SP depending on whether the middle SP is
      // in the central or forward region
      nTopSeedConf = rM > seedConfRange.rMaxSeedConf
                         ? seedConfRange.nTopForLargeR
                         : seedConfRange.nTopForSmallR;
    }

    state.compatTopSP.clear();

//    std::cout << "======= MIDDLE " << rM << " =======" << std::endl;

//    std::cout << "TEST TOP" << std::endl;

    for (auto topSP : topSPs) {
      float rT = topSP->radius();
      float deltaR = rT - rM;

//      std::cout << "rT " << rT << " " << deltaR << std::endl;

      // if r-distance is too small, try next SP in bin
      if (deltaR < m_config.deltaRMinTopSP) {
        continue;
      }
      // if r-distance is too big, try next SP in bin
      if (deltaR > m_config.deltaRMaxTopSP) {
        // if SPs are sorted in radius we break out of the loop
        if (m_config.forceRadialSorting) {
//					std::cout << "BREAK" << std::endl;
          break;
        }
        continue;
      }
      float deltaZ = topSP->z() - zM;
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = deltaZ / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        continue;
      }
      if (std::abs(deltaZ) > m_config.deltaZMax) {
        continue;
      }
      // cut on the max curvature between top SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      if (m_config.interactionPointCut) {
        float xVal = (topSP->x() - spM->x()) * (spM->x() / rM) +
                     (topSP->y() - spM->y()) * (spM->y() / rM);
        float yVal = (topSP->y() - spM->y()) * (spM->x() / rM) -
                     (topSP->x() - spM->x()) * (spM->y() / rM);
        if (std::abs(rM * yVal) > m_config.impactMax * xVal) {
          // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
          // circle into straight lines in the u/v plane the line equation can
          // be described in terms of aCoef and bCoef, where v = aCoef * u +
          // bCoef
          float uT = xVal / (xVal * xVal + yVal * yVal);
          float vT = yVal / (xVal * xVal + yVal * yVal);
          // in the rotated frame the interaction point is positioned at x = -rM
          // and y ~= impactParam
          float uIP = -1. / rM;
          float vIP = m_config.impactMax / (rM * rM);
          if (yVal > 0.)
            vIP = -vIP;
          // we can obtain aCoef as the slope dv/du of the linear function,
          // estimated using du and dv between the two SP bCoef is obtained by
          // inserting aCoef into the linear equation
          float aCoef = (vT - vIP) / (uT - uIP);
          float bCoef = vIP - aCoef * uIP;
          // the distance of the straight line from the origin (radius of the
          // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
          // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
          if ((bCoef * bCoef) >
              (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
            continue;
          }
        }
      }
      state.compatTopSP.push_back(topSP);
    }
    if (state.compatTopSP.empty()) {
      continue;
    }
    // apply cut on the number of top SP if seedConfirmation is true
    if (m_config.seedConfirmation == true &&
        state.compatTopSP.size() < nTopSeedConf) {
      continue;
    }

    state.compatBottomSP.clear();

//    std::cout << "TEST BOT" << std::endl;
	
    for (auto bottomSP : bottomSPs) {
      float rB = bottomSP->radius();
      float deltaR = rM - rB;

//      std::cout << "rB "  << rB << " " << deltaR << std::endl;

      // this condition is the opposite of the condition for top SP
      if (deltaR > m_config.deltaRMaxBottomSP) {
        continue;
      }
      if (deltaR < m_config.deltaRMinBottomSP) {
        // if SPs are sorted in radius we break out of the loop
        if (m_config.forceRadialSorting) {
//					std::cout << "BREAK" << std::endl;
          break;
        }
        continue;
      }
      float deltaZ = zM - bottomSP->z();
      // ratio Z/R (forward angle) of space point duplet
      float cotTheta = deltaZ / deltaR;
      if (std::fabs(cotTheta) > m_config.cotThetaMax) {
        continue;
      }
      // check if duplet origin on z axis within collision region
      float zOrigin = zM - rM * cotTheta;
      if (zOrigin < m_config.collisionRegionMin ||
          zOrigin > m_config.collisionRegionMax) {
        continue;
      }
      if (std::abs(deltaZ) > m_config.deltaZMax) {
        continue;
      }
      // cut on the max curvature between bottom SP and interaction point
      // first transform the space point coordinates into a frame such that the
      // central space point SPm is in the origin of the frame and the x axis
      // points away from the interaction point in addition to a translation
      // transformation we also perform a rotation in order to keep the
      // curvature of the circle tangent to the x axis
      if (m_config.interactionPointCut) {
        float xVal = (bottomSP->x() - spM->x()) * (spM->x() / rM) +
                     (bottomSP->y() - spM->y()) * (spM->y() / rM);
        float yVal = (bottomSP->y() - spM->y()) * (spM->x() / rM) -
                     (bottomSP->x() - spM->x()) * (spM->y() / rM);
        if (std::abs(rM * yVal) > -m_config.impactMax * xVal) {
          // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
          // circle into straight lines in the u/v plane the line equation can
          // be
          // described in terms of aCoef and bCoef, where v = aCoef * u + bCoef
          float uB = xVal / (xVal * xVal + yVal * yVal);
          float vB = yVal / (xVal * xVal + yVal * yVal);
          // in the rotated frame the interaction point is positioned at x = -rM
          // and y ~= impactParam
          float uIP = -1. / rM;
          float vIP = m_config.impactMax / (rM * rM);
          if (yVal < 0.)
            vIP = -vIP;
          // we can obtain aCoef as the slope dv/du of the linear function,
          // estimated using du and dv between the two SP bCoef is obtained by
          // inserting aCoef into the linear equation
          float aCoef = (vB - vIP) / (uB - uIP);
          float bCoef = vIP - aCoef * uIP;
          // the distance of the straight line from the origin (radius of the
          // circle) is related to aCoef and bCoef by d^2 = bCoef^2 / (1 +
          // aCoef^2) = 1 / (radius^2) and we can apply the cut on the curvature
          if ((bCoef * bCoef) >
              (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
            continue;
          }
        }
      }
      state.compatBottomSP.push_back(bottomSP);
    }
    // no bottom SP found -> try next spM
    if (state.compatBottomSP.empty()) {
      continue;
    }

    state.linCircleBottom.clear();
    state.linCircleTop.clear();

    transformCoordinates(state.compatBottomSP, *spM, true,
                         state.linCircleBottom);
    transformCoordinates(state.compatTopSP, *spM, false, state.linCircleTop);

    state.topSpVec.clear();
    state.curvatures.clear();
    state.impactParameters.clear();
    state.seedsPerSpM.clear();

    size_t numBotSP = state.compatBottomSP.size();
    size_t numTopSP = state.compatTopSP.size();

    int numQualitySeeds = 0;
    int numSeeds = 0;

    size_t t0 = 0;

    for (size_t b = 0; b < numBotSP; b++) {
      auto lb = state.linCircleBottom[b];
      float Zob = lb.Zo;
      float cotThetaB = lb.cotTheta;
      float Vb = lb.V;
      float Ub = lb.U;
      float ErB = lb.Er;
      float iDeltaRB = lb.iDeltaR;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1. + cotThetaB * cotThetaB);
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
      // scattering
      // but to avoid trig functions we approximate cot by scaling by
      // 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at
      // eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *=
          m_config.sigmaScattering * m_config.sigmaScattering;

      float sinTheta = 1 / std::sqrt(iSinTheta2);
      float cosTheta = cotThetaB * sinTheta;

      // clear all vectors used in each inner for loop
      state.topSpVec.clear();
      state.curvatures.clear();
      state.impactParameters.clear();
      for (size_t t = t0; t < numTopSP; t++) {
        auto lt = state.linCircleTop[t];

        float cotThetaT = lt.cotTheta;
        float rMxy = 0.;
        float ub = 0.;
        float vb = 0.;
        float ut = 0.;
        float vt = 0.;

        if (m_config.useDetailedDoubleMeasurementInfo) {
          // protects against division by 0
          float dU = lt.U - Ub;
          if (dU == 0.) {
            continue;
          }
          // A and B are evaluated as a function of the circumference parameters
          // x_0 and y_0
          float A0 = (lt.V - Vb) / dU;

          // coordinate transformation and checks for middle spacepoint
          // x and y terms for the rotation from UV to XY plane
          float rotationTermsUVtoXY[2] = {spM->x() * sinTheta / spM->radius(),
                                          spM->y() * sinTheta / spM->radius()};
          // position of Middle SP converted from UV to XY assuming cotTheta
          // evaluated from the Bottom and Middle SPs double
          double positionMiddle[3] = {
              rotationTermsUVtoXY[0] - rotationTermsUVtoXY[1] * A0,
              rotationTermsUVtoXY[0] * A0 + rotationTermsUVtoXY[1],
              cosTheta * std::sqrt(1 + A0 * A0)};

          double rMTransf[3];
          if (!xyzCoordinateCheck(m_config, spM, positionMiddle,
                                  m_config.toleranceParam, rMTransf)) {
            continue;
          }

          // coordinate transformation and checks for bottom spacepoint
          float B0 = 2. * (Vb - A0 * Ub);
          float Cb = 1. - B0 * lb.y;
          float Sb = A0 + B0 * lb.x;
          double positionBottom[3] = {
              rotationTermsUVtoXY[0] * Cb - rotationTermsUVtoXY[1] * Sb,
              rotationTermsUVtoXY[0] * Sb + rotationTermsUVtoXY[1] * Cb,
              cosTheta * std::sqrt(1 + A0 * A0)};

          auto spB = state.compatBottomSP[b];
          double rBTransf[3];
          if (!xyzCoordinateCheck(m_config, spB, positionBottom,
                                  m_config.toleranceParam, rBTransf)) {
            continue;
          }

          // coordinate transformation and checks for top spacepoint
          float Ct = 1. - B0 * lt.y;
          float St = A0 + B0 * lt.x;
          double positionTop[3] = {
              rotationTermsUVtoXY[0] * Ct - rotationTermsUVtoXY[1] * St,
              rotationTermsUVtoXY[0] * St + rotationTermsUVtoXY[1] * Ct,
              cosTheta * std::sqrt(1 + A0 * A0)};

          auto spT = state.compatTopSP[t];
          double rTTransf[3];
          if (!xyzCoordinateCheck(m_config, spT, positionTop,
                                  m_config.toleranceParam, rTTransf)) {
            continue;
          }

          // bottom and top coordinates in the spM reference frame
          float xB = rBTransf[0] - rMTransf[0];
          float yB = rBTransf[1] - rMTransf[1];
          float zB = rBTransf[2] - rMTransf[2];
          float xT = rTTransf[0] - rMTransf[0];
          float yT = rTTransf[1] - rMTransf[1];
          float zT = rTTransf[2] - rMTransf[2];

          float iDeltaRB2 = 1. / (xB * xB + yB * yB);
          float iDeltaRT2 = 1. / (xT * xT + yT * yT);

          cotThetaB = -zB * std::sqrt(iDeltaRB2);
          cotThetaT = zT * std::sqrt(iDeltaRT2);

          rMxy =
              std::sqrt(rMTransf[0] * rMTransf[0] + rMTransf[1] * rMTransf[1]);
          float Ax = rMTransf[0] / rMxy;
          float Ay = rMTransf[1] / rMxy;

          ub = (xB * Ax + yB * Ay) * iDeltaRB2;
          vb = (yB * Ax - xB * Ay) * iDeltaRB2;
          ut = (xT * Ax + yT * Ay) * iDeltaRT2;
          vt = (yT * Ax - xT * Ay) * iDeltaRT2;
        }

        // use geometric average
        float cotThetaAvg2 = cotThetaB * cotThetaT;
        if (m_config.arithmeticAverageCotTheta) {
          // use arithmetic average
          cotThetaAvg2 = std::pow((cotThetaB + cotThetaT) / 2, 2);
        }

        // add errors of spB-spM and spM-spT pairs and add the correlation term
        // for errors on spM
        float error2 = lt.Er + ErB +
                       2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB *
                           lt.iDeltaR;

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
        if (deltaCotTheta2 > (error2 + scatteringInRegion2)) {
          // additional cut to skip top SPs when producing triplets
          if (m_config.skipPreviousTopSP) {
            // break if cotTheta from bottom SP < cotTheta from top SP because
            // the SP are sorted by cotTheta
            if (cotThetaB - cotThetaT < 0) {
              break;
            }
            t0 = t + 1;
          }
          continue;
        }

        float dU;
        float A;
        float S2;
        float B;
        float B2;

        if (m_config.useDetailedDoubleMeasurementInfo) {
          dU = ut - ub;
          // protects against division by 0
          if (dU == 0.) {
            continue;
          }
          A = (vt - vb) / dU;
          S2 = 1. + A * A;
          B = vb - A * ub;
          B2 = B * B;
        } else {
          dU = lt.U - Ub;
          // protects against division by 0
          if (dU == 0.) {
            continue;
          }
          // A and B are evaluated as a function of the circumference parameters
          // x_0 and y_0
          A = (lt.V - Vb) / dU;
          S2 = 1. + A * A;
          B = Vb - A * Ub;
          B2 = B * B;
        }

        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        if (S2 < B2 * m_config.minHelixDiameter2) {
          continue;
        }

        // refinement of the cut on the compatibility between the r-z slope of
        // the two seed segments using a scattering term scaled by the actual
        // measured pT (p2scatterSigma)
        float iHelixDiameter2 = B2 / S2;
        // calculate scattering for p(T) calculated from seed curvature
        float pT2scatterSigma = iHelixDiameter2 * m_config.sigmapT2perRadius;
        // if pT > maxPtScattering, calculate allowed scattering angle using
        // maxPtScattering instead of pt.
        float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
        if (pT > m_config.maxPtScattering) {
          float pTscatterSigma =
              (m_config.highland / m_config.maxPtScattering) *
              m_config.sigmaScattering;
          pT2scatterSigma = pTscatterSigma * pTscatterSigma;
        }
        // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
        // from rad to deltaCotTheta
        float p2scatterSigma = pT2scatterSigma * iSinTheta2;
        // if deltaTheta larger than allowed scattering for calculated pT, skip
        if (deltaCotTheta2 > (error2 + p2scatterSigma)) {
          if (m_config.skipPreviousTopSP) {
            if (cotThetaB - cotThetaT < 0) {
              break;
            }
            t0 = t;
          }
          continue;
        }
        // A and B allow calculation of impact params in U/V plane with linear
        // function
        // (in contrast to having to solve a quadratic function in x/y plane)
        float Im;
        if (m_config.useDetailedDoubleMeasurementInfo == false) {
          Im = std::abs((A - B * rM) * rM);
        } else {
          Im = std::abs((A - B * rMxy) * rMxy);
        }

        if (Im <= m_config.impactMax) {
          state.topSpVec.push_back(state.compatTopSP[t]);
          // inverse diameter is signed depending if the curvature is
          // positive/negative in phi
          state.curvatures.push_back(B / std::sqrt(S2));
          state.impactParameters.push_back(Im);

          // evaluate eta and pT of the seed
          float cotThetaAvg = std::sqrt(cotThetaAvg2);
          float theta = std::atan(1. / cotThetaAvg);
          float eta = -std::log(std::tan(0.5 * theta));
          state.etaVec.push_back(eta);
          state.ptVec.push_back(pT);
        }
      }
      if (!state.topSpVec.empty()) {
        m_config.seedFilter->filterSeeds_2SpFixed(
            *state.compatBottomSP[b], *spM, state.topSpVec, state.curvatures,
            state.impactParameters, Zob, numQualitySeeds, numSeeds,
            state.seedsPerSpM);
      }
    }
    m_config.seedFilter->filterSeeds_1SpFixed(state.seedsPerSpM,
                                              numQualitySeeds, outIt);
  }
}

template <typename external_spacepoint_t, typename platform_t>
template <typename sp_range_t>
std::vector<Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  State state;
  Extent extent;
  std::vector<Seed<external_spacepoint_t>> ret;

  createSeedsForGroup(state, std::back_inserter(ret), bottomSPs, middleSPs,
                      topSPs, extent);

  return ret;
}
}  // namespace Acts
