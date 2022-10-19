// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"

#include <cmath>
#include <functional>
#include <numeric>
#include <type_traits>

namespace Acts {
template <typename external_spacepoint_t>
auto SeedFinderOrthogonal<external_spacepoint_t>::validTupleOrthoRangeLH(
    const internal_sp_t &low) const -> typename tree_t::range_t {
  float colMin = m_config.collisionRegionMin;
  float colMax = m_config.collisionRegionMax;
  float pL = low.phi();
  float rL = low.radius();
  float zL = low.z();

  typename tree_t::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(m_config.phiMin);
  res[DimPhi].shrinkMax(m_config.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(m_config.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(m_config.zMin);
  res[DimZ].shrinkMax(m_config.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r - r_L ≤ Δr_max, as defined
   * by the seeding configuration and the given lower spacepoint.
   */
  res[DimR].shrinkMin(rL + m_config.deltaRMinTopSP);
  res[DimR].shrinkMax(rL + m_config.deltaRMaxTopSP);

  /*
   * Cut: Now that we have constrained r, we can use that new r range to
   * further constrain z.
   */
  float zMax = (res[DimR].max() / rL) * (zL - colMin) + colMin;
  float zMin = colMax - (res[DimR].max() / rL) * (colMax - zL);

  /*
   * This cut only works if z_low is outside the collision region for z.
   */
  if (zL > colMin) {
    res[DimZ].shrinkMax(zMax);
  } else if (zL < colMax) {
    res[DimZ].shrinkMin(zMin);
  }

  /*
   * Cut: Shrink the z-range using the maximum cot(θ).
   */
  res[DimZ].shrinkMin(zL - m_config.cotThetaMax * (res[DimR].max() - rL));
  res[DimZ].shrinkMax(zL + m_config.cotThetaMax * (res[DimR].max() - rL));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_L ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pL - m_config.deltaPhiMax);
  res[DimPhi].shrinkMax(pL + m_config.deltaPhiMax);

  return res;
}

template <typename external_spacepoint_t>
auto SeedFinderOrthogonal<external_spacepoint_t>::validTupleOrthoRangeHL(
    const internal_sp_t &high) const -> typename tree_t::range_t {
  float pM = high.phi();
  float rM = high.radius();

  typename tree_t::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(m_config.phiMin);
  res[DimPhi].shrinkMax(m_config.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(m_config.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(m_config.zMin);
  res[DimZ].shrinkMax(m_config.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r_H - r ≤ Δr_max, as defined
   * by the seeding configuration and the given higher spacepoint.
   */
  res[DimR].shrinkMin(rM - m_config.deltaRMaxBottomSP);
  res[DimR].shrinkMax(rM - m_config.deltaRMinBottomSP);

  /*
   * Cut: Now that we have constrained r, we can use that new r range to
   * further constrain z.
   */
  float fracR = res[DimR].min() / rM;

  float zMin = (high.z() - m_config.collisionRegionMin) * fracR +
               m_config.collisionRegionMin;
  float zMax = (high.z() - m_config.collisionRegionMax) * fracR +
               m_config.collisionRegionMax;

  res[DimZ].shrinkMin(std::min(zMin, high.z()));
  res[DimZ].shrinkMax(std::max(zMax, high.z()));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_H ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pM - m_config.deltaPhiMax);
  res[DimPhi].shrinkMax(pM + m_config.deltaPhiMax);

  return res;
}

template <typename external_spacepoint_t>
bool SeedFinderOrthogonal<external_spacepoint_t>::validTuple(
    const internal_sp_t &low, const internal_sp_t &high) const {
  float rL = low.radius();
  float rH = high.radius();

  float zL = low.z();
  float zH = high.z();

  float deltaR = rH - rL;

  /*
   * Cut: Ensure that the forward angle (z / r) lies within reasonable bounds,
   * which is to say the absolute value must be smaller than the max cot(θ).
   */
  float deltaZ = (zH - zL);
  float cotTheta = deltaZ / deltaR;
  if (std::fabs(cotTheta) > m_config.cotThetaMax) {
    return false;
  }

  /*
   * Cut: Ensure that the origin of the duplet (the intersection of the line
   * between them with the z axis) lies within the collision region.
   */
  float zOrigin = zL - rL * cotTheta;
  if (zOrigin < m_config.collisionRegionMin ||
      zOrigin > m_config.collisionRegionMax) {
    return false;
  }
  if (std::abs(deltaZ) > m_config.deltaZMax) {
    return false;
  }

  // cut on the max curvature calculated from dublets
  // first transform the space point coordinates into a frame such that the
  // central space point SPm is in the origin of the frame and the x axis
  // points away from the interaction point in addition to a translation
  // transformation we also perform a rotation in order to keep the
  // curvature of the circle tangent to the x axis
  if (m_config.interactionPointCut) {
    float xVal = (high.x() - low.x()) * (low.x() / rL) +
                 (high.y() - low.y()) * (low.y() / rL);
    float yVal = (high.y() - low.y()) * (low.x() / rL) -
                 (high.x() - low.x()) * (low.y() / rL);
    if (std::abs(rL * yVal) > m_config.impactMax * xVal) {
      // conformal transformation u=x/(x²+y²) v=y/(x²+y²) transform the
      // circle into straight lines in the u/v plane the line equation can
      // be described in terms of aCoef and bCoef, where v = aCoef * u +
      // bCoef
      float uT = xVal / (xVal * xVal + yVal * yVal);
      float vT = yVal / (xVal * xVal + yVal * yVal);
      // in the rotated frame the interaction point is positioned at x = -rM
      // and y ~= impactParam
      float uIP = -1. / rL;
      float vIP = m_config.impactMax / (rL * rL);
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
      if ((bCoef * bCoef) > (1 + aCoef * aCoef) / m_config.minHelixDiameter2) {
        return false;
      }
    }
  }

  return true;
}

template <typename external_spacepoint_t>
SeedFinderOrthogonal<external_spacepoint_t>::SeedFinderOrthogonal(
    SeedFinderOrthogonalConfig<external_spacepoint_t> config)
    : m_config(config.toInternalUnits()) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(config.radLengthPerSeed));
  float maxScatteringAngle = config.highland / config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(config.minPt * 2 / config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(config.highland / config.pTPerHelixRadius, 2);
  m_config.sigmapT2perRadius =
      config.pT2perRadius * std::pow(2 * config.sigmaScattering, 2);
}

template <typename external_spacepoint_t>
template <typename output_container_t>
void SeedFinderOrthogonal<external_spacepoint_t>::filterCandidates(
    internal_sp_t &middle, std::vector<internal_sp_t *> &bottom,
    std::vector<internal_sp_t *> &top, SeedFilterState seedFilterState,
    output_container_t &cont) const {
  float rM = middle.radius();
  float varianceRM = middle.varianceR();
  float varianceZM = middle.varianceZ();

  std::vector<internal_sp_t *> top_valid;
  std::vector<float> curvatures;
  std::vector<float> impactParameters;

  // contains parameters required to calculate circle with linear equation
  // ...for bottom-middle
  std::vector<LinCircle> linCircleBottom;
  linCircleBottom.reserve(bottom.size());
  // ...for middle-top
  std::vector<LinCircle> linCircleTop;
  linCircleTop.reserve(top.size());

  transformCoordinates(bottom, middle, true, linCircleBottom);
  transformCoordinates(top, middle, false, linCircleTop);

  std::vector<float> tanLM;
  std::vector<float> tanMT;

  tanLM.reserve(bottom.size());
  tanMT.reserve(top.size());

  size_t numBotSP = bottom.size();
  size_t numTopSP = top.size();

  for (size_t b = 0; b < numBotSP; b++) {
    tanLM.push_back(std::atan2(middle.radius() - bottom[b]->radius(),
                               middle.z() - bottom[b]->z()));
  }

  for (size_t t = 0; t < numTopSP; t++) {
    tanMT.push_back(std::atan2(top[t]->radius() - middle.radius(),
                               top[t]->z() - middle.z()));
  }

  for (size_t b = 0; b < numBotSP; b++) {
    auto lb = linCircleBottom[b];
    seedFilterState.zOrigin = lb.Zo;
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
    scatteringInRegion2 *= m_config.sigmaScattering * m_config.sigmaScattering;

    // clear all vectors used in each inner for loop
    top_valid.clear();
    curvatures.clear();
    impactParameters.clear();
    for (size_t t = 0; t < numTopSP; t++) {
      auto lt = linCircleTop[t];

      if (std::abs(tanLM[b] - tanMT[t]) > 0.005) {
        continue;
      }

      // add errors of spB-spM and spM-spT pairs and add the correlation term
      // for errors on spM
      float error2 = lt.Er + ErB +
                     2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
                         iDeltaRB * lt.iDeltaR;

      float deltaCotTheta = cotThetaB - lt.cotTheta;
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
        continue;
      }

      float dU = lt.U - Ub;

      // A and B are evaluated as a function of the circumference parameters
      // x_0 and y_0
      float A = (lt.V - Vb) / dU;
      float S2 = 1. + A * A;
      float B = Vb - A * Ub;
      float B2 = B * B;
      // sqrt(S2)/B = 2 * helixradius
      // calculated radius must not be smaller than minimum radius
      if (S2 < B2 * m_config.minHelixDiameter2) {
        continue;
      }
      // 1/helixradius: (B/sqrt(S2))*2 (we leave everything squared)
      float iHelixDiameter2 = B2 / S2;
      // calculate scattering for p(T) calculated from seed curvature
      float pT2scatter = 4 * iHelixDiameter2 * m_config.pT2perRadius;
      // if pT > maxPtScattering, calculate allowed scattering angle using
      // maxPtScattering instead of pt.
      float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
      if (pT > m_config.maxPtScattering) {
        float pTscatter = m_config.highland / m_config.maxPtScattering;
        pT2scatter = pTscatter * pTscatter;
      }
      // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
      // from rad to deltaCotTheta
      // if deltaTheta larger than allowed scattering for calculated pT, skip
      if (deltaCotTheta2 >
          (error2 + (pT2scatter * iSinTheta2 * m_config.sigmaScattering *
                     m_config.sigmaScattering))) {
        continue;
      }

      // A and B allow calculation of impact params in U/V plane with linear
      // function
      // (in contrast to having to solve a quadratic function in x/y plane)
      float Im = std::abs((A - B * rM) * rM);

      if (Im <= m_config.impactMax) {
        top_valid.push_back(top[t]);
        // inverse diameter is signed depending if the curvature is
        // positive/negative in phi
        curvatures.push_back(B / std::sqrt(S2));
        impactParameters.push_back(Im);
      }
    }
    if (!top_valid.empty()) {
      m_config.seedFilter->filterSeeds_2SpFixed(*bottom[b], middle, top_valid,
                                                curvatures, impactParameters,
                                                seedFilterState, cont);
    }
  }
}

template <typename external_spacepoint_t>
template <typename output_container_t>
void SeedFinderOrthogonal<external_spacepoint_t>::processFromMiddleSP(
    const tree_t &tree, output_container_t &out_cont,
    const typename tree_t::pair_t &middle_p) const {
  using range_t = typename tree_t::range_t;
  internal_sp_t &middle = *middle_p.second;

  /*
   * Prepare four output vectors for seed candidates:
   *
   * bottom_lh_v denotes the candidates bottom seed points, assuming that the
   * track has monotonically _increasing_ z position. bottom_hl_v denotes the
   * candidate bottom points assuming that the track has monotonically
   * _decreaing_ z position. top_lh_v are the candidate top points for an
   * increasing z track, and top_hl_v are the candidate top points for a
   * decreasing z track.
   */
  std::vector<internal_sp_t *> bottom_lh_v, bottom_hl_v, top_lh_v, top_hl_v;

  /*
   * Cut: Ensure that the middle spacepoint lies within a valid r-region for
   * middle points.
   */
  if (middle.radius() > m_config.rMaxMiddle ||
      middle.radius() < m_config.rMinMiddle) {
    return;
  }

  /*
   * Calculate the search ranges for bottom and top candidates for this middle
   * space point.
   */
  range_t bottom_r = validTupleOrthoRangeHL(middle);
  range_t top_r = validTupleOrthoRangeLH(middle);

  /*
   * Calculate the value of cot(θ) for this middle spacepoint.
   */
  float myCotTheta =
      std::max(std::abs(middle.z() / middle.radius()), m_config.cotThetaMax);

  /*
   * Calculate the maximum Δr, given that we have already constrained our
   * search space.
   */
  float deltaRMaxTop = top_r[DimR].max() - middle.radius();
  float deltaRMaxBottom = middle.radius() - bottom_r[DimR].min();

  /*
   * Create the search range for the bottom spacepoint assuming a monotonically
   * increasing z track, by calculating the minimum z value from the cot(θ),
   * and by setting the maximum to the z position of the middle spacepoint - if
   * the z position is higher than the middle point, then it would be a
   * decreasing z track!
   */
  range_t bottom_lh_r = bottom_r;
  bottom_lh_r[DimZ].shrink(middle.z() - myCotTheta * deltaRMaxBottom,
                           middle.z());

  /*
   * Calculate the search ranges for the other four sets of points in a similar
   * fashion.
   */
  range_t top_lh_r = top_r;
  top_lh_r[DimZ].shrink(middle.z(), middle.z() + myCotTheta * deltaRMaxTop);

  range_t bottom_hl_r = bottom_r;
  bottom_hl_r[DimZ].shrink(middle.z(),
                           middle.z() + myCotTheta * deltaRMaxBottom);
  range_t top_hl_r = top_r;
  top_hl_r[DimZ].shrink(middle.z() - myCotTheta * deltaRMaxTop, middle.z());

  /*
   * Make sure the candidate vectors are clear, in case we've used them before.
   */
  bottom_lh_v.clear();
  bottom_hl_v.clear();
  top_lh_v.clear();
  top_hl_v.clear();

  /*
   * Now, we will actually search for the spaces. Remembering that we combine
   * bottom and top candidates for increasing and decreasing tracks separately,
   * we will first check whether both the search ranges for increasing tracks
   * are not degenerate - if they are, we will never find any seeds and we do
   * not need to bother doing the search.
   */
  if (!bottom_lh_r.degenerate() && !top_lh_r.degenerate()) {
    /*
     * Search the trees for points that lie in the given search range.
     */
    tree.rangeSearchMapDiscard(
        bottom_lh_r,
        [this, &middle, &bottom_lh_v](const typename tree_t::coordinate_t &,
                                      const typename tree_t::value_t &bottom) {
          if (validTuple(*bottom, middle)) {
            bottom_lh_v.push_back(bottom);
          }
        });
  }

  /*
   * Perform the same search for candidate bottom spacepoints, but for
   * monotonically decreasing z tracks.
   */
  if (!bottom_hl_r.degenerate() && !top_hl_r.degenerate()) {
    tree.rangeSearchMapDiscard(
        bottom_hl_r,
        [this, &middle, &bottom_hl_v](const typename tree_t::coordinate_t &,
                                      const typename tree_t::value_t &bottom) {
          if (validTuple(middle, *bottom)) {
            bottom_hl_v.push_back(bottom);
          }
        });
  }

  /*
   * Next, we perform a search for top candidates in increasing z tracks, which
   * only makes sense if we found any bottom candidates.
   */
  if (!bottom_lh_v.empty()) {
    tree.rangeSearchMapDiscard(
        top_lh_r,
        [this, &middle, &top_lh_v](const typename tree_t::coordinate_t &,
                                   const typename tree_t::value_t &top) {
          if (validTuple(*top, middle)) {
            top_lh_v.push_back(top);
          }
        });
  }

  /*
   * And repeat for the top spacepoints for decreasing z tracks!
   */
  if (!bottom_hl_v.empty()) {
    tree.rangeSearchMapDiscard(
        top_hl_r,
        [this, &middle, &top_hl_v](const typename tree_t::coordinate_t &,
                                   const typename tree_t::value_t &top) {
          if (validTuple(middle, *top)) {
            top_hl_v.push_back(top);
          }
        });
  }

  /*
   * Create a vector to contain protoseeds.
   */
  std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
      protoseeds;

  // TODO: add seed confirmation
  SeedFilterState seedFilterState;

  /*
   * If we have candidates for increasing z tracks, we try to combine them.
   */
  if (!bottom_lh_v.empty() && !top_lh_v.empty()) {
    filterCandidates(middle, bottom_lh_v, top_lh_v, seedFilterState,
                     protoseeds);
  }

  /*
   * Try to combine candidates for decreasing z tracks.
   */
  if (!bottom_hl_v.empty() && !top_hl_v.empty()) {
    filterCandidates(middle, bottom_hl_v, top_hl_v, seedFilterState,
                     protoseeds);
  }

  /*
   * Run a seed filter, just like in other seeding algorithms.
   */
  m_config.seedFilter->filterSeeds_1SpFixed(protoseeds,
                                            seedFilterState.numQualitySeeds,
                                            std::back_inserter(out_cont));
}

template <typename external_spacepoint_t>
auto SeedFinderOrthogonal<external_spacepoint_t>::createTree(
    const std::vector<internal_sp_t *> &spacePoints) const -> tree_t {
  std::vector<typename tree_t::pair_t> points;

  /*
   * For every input point, we create a coordinate-pointer pair, which we then
   * linearly pass to the k-d tree constructor. That constructor will take care
   * of sorting the pairs and splitting the space.
   */
  for (internal_sp_t *sp : spacePoints) {
    typename tree_t::coordinate_t point;

    point[DimPhi] = sp->phi();
    point[DimR] = sp->radius();
    point[DimZ] = sp->z();

    points.emplace_back(point, sp);
  }

  return tree_t(std::move(points));
}

template <typename external_spacepoint_t>
template <typename input_container_t, typename output_container_t>
void SeedFinderOrthogonal<external_spacepoint_t>::createSeeds(
    const input_container_t &spacePoints, output_container_t &out_cont) const {
  /*
   * The template parameters we accept are a little too generic, so we want to
   * run some basic checks to make sure the containers have the correct value
   * types.
   */
  static_assert(std::is_same_v<typename output_container_t::value_type,
                               Seed<external_spacepoint_t>>,
                "Output iterator container type must accept seeds.");
  static_assert(std::is_same_v<typename input_container_t::value_type,
                               const external_spacepoint_t *>,
                "Input container must contain external spacepoints.");

  /*
   * Sadly, for the time being, we will need to construct our internal space
   * points on the heap. This adds some additional overhead and work. Here we
   * take each external spacepoint, allocate a corresponding internal space
   * point, and save it in a vector.
   */
  std::vector<internal_sp_t *> internalSpacePoints;
  for (const external_spacepoint_t *p : spacePoints) {
    internalSpacePoints.push_back(new InternalSpacePoint<external_spacepoint_t>(
        *p, {p->x(), p->y(), p->z()}, {0.0, 0.0},
        {p->varianceR(), p->varianceZ()}));
  }

  /*
   * Construct the k-d tree from these points. Note that this not consume or
   * take ownership of the points.
   */
  tree_t tree = createTree(internalSpacePoints);

  /*
   * Run the seeding algorithm by iterating over all the points in the tree and
   * seeing what happens if we take them to be our middle spacepoint.
   */
  for (const typename tree_t::pair_t &middle_p : tree) {
    processFromMiddleSP(tree, out_cont, middle_p);
  }

  /*
   * Don't forget to get rid of all the spacepoints we just allocated!
   */
  for (const internal_sp_t *p : internalSpacePoints) {
    delete p;
  }
}

template <typename external_spacepoint_t>
template <typename input_container_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderOrthogonal<external_spacepoint_t>::createSeeds(
    const input_container_t &spacePoints) const {
  std::vector<seed_t> r;

  createSeeds(spacePoints, r);

  return r;
}

}  // namespace Acts
