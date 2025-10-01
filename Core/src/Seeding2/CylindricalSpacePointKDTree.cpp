// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/CylindricalSpacePointKDTree.hpp"

namespace Acts::Experimental {

CylindricalSpacePointKDTree::CylindricalSpacePointKDTree(
    Tree tree, std::unique_ptr<const Logger> logger)
    : m_tree(std::move(tree)), m_logger(std::move(logger)) {}

CylindricalSpacePointKDTree::Tree::range_t
CylindricalSpacePointKDTree::validTupleOrthoRangeLH(
    const Options &options, const ConstSpacePointProxy2 &low) const {
  float colMin = options.collisionRegionMin;
  float colMax = options.collisionRegionMax;
  float pL = low.phi();
  float rL = low.zr()[1];
  float zL = low.zr()[0];

  Tree::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(options.phiMin);
  res[DimPhi].shrinkMax(options.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(options.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(options.zMin);
  res[DimZ].shrinkMax(options.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r - r_L ≤ Δr_max, as defined
   * by the seeding configuration and the given lower spacepoint.
   */
  res[DimR].shrinkMin(rL + options.deltaRMin);
  res[DimR].shrinkMax(rL + options.deltaRMax);

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
  res[DimZ].shrinkMin(zL - options.cotThetaMax * (res[DimR].max() - rL));
  res[DimZ].shrinkMax(zL + options.cotThetaMax * (res[DimR].max() - rL));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_L ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pL - options.deltaPhiMax);
  res[DimPhi].shrinkMax(pL + options.deltaPhiMax);

  // Cut: Ensure that z-distance between SPs is within max and min values.
  res[DimZ].shrinkMin(zL - options.deltaZMax);
  res[DimZ].shrinkMax(zL + options.deltaZMax);

  return res;
}

CylindricalSpacePointKDTree::Tree::range_t
CylindricalSpacePointKDTree::validTupleOrthoRangeHL(
    const Options &options, const ConstSpacePointProxy2 &high) const {
  float pM = high.phi();
  float rM = high.zr()[1];
  float zM = high.zr()[0];

  Tree::range_t res;

  /*
   * Cut: Ensure that we search only in φ_min ≤ φ ≤ φ_max, as defined by the
   * seeding configuration.
   */
  res[DimPhi].shrinkMin(options.phiMin);
  res[DimPhi].shrinkMax(options.phiMax);

  /*
   * Cut: Ensure that we search only in r ≤ r_max, as defined by the seeding
   * configuration.
   */
  res[DimR].shrinkMax(options.rMax);

  /*
   * Cut: Ensure that we search only in z_min ≤ z ≤ z_max, as defined by the
   * seeding configuration.
   */
  res[DimZ].shrinkMin(options.zMin);
  res[DimZ].shrinkMax(options.zMax);

  /*
   * Cut: Ensure that we search only in Δr_min ≤ r_H - r ≤ Δr_max, as defined
   * by the seeding configuration and the given higher spacepoint.
   */
  res[DimR].shrinkMin(rM - options.deltaRMax);
  res[DimR].shrinkMax(rM - options.deltaRMin);

  /*
   * Cut: Now that we have constrained r, we can use that new r range to
   * further constrain z.
   */
  float fracR = res[DimR].min() / rM;

  float zMin =
      (zM - options.collisionRegionMin) * fracR + options.collisionRegionMin;
  float zMax =
      (zM - options.collisionRegionMax) * fracR + options.collisionRegionMax;

  res[DimZ].shrinkMin(std::min(zMin, zM));
  res[DimZ].shrinkMax(std::max(zMax, zM));

  /*
   * Cut: Shrink the φ range, such that Δφ_min ≤ φ - φ_H ≤ Δφ_max
   */
  res[DimPhi].shrinkMin(pM - options.deltaPhiMax);
  res[DimPhi].shrinkMax(pM + options.deltaPhiMax);

  // Cut: Ensure that z-distance between SPs is within max and min values.
  res[DimZ].shrinkMin(zM - options.deltaZMax);
  res[DimZ].shrinkMax(zM + options.deltaZMax);

  return res;
}

void CylindricalSpacePointKDTree::validTuples(const Options &lhOptions,
                                              const Options &hlOptions,
                                              const ConstSpacePointProxy2 &spM,
                                              std::size_t nTopSeedConf,
                                              Candidates &candidates) const {
  using range_t = Tree::range_t;

  /*
   * Calculate the search ranges for bottom and top candidates for this middle
   * space point.
   */
  range_t bottom_r = validTupleOrthoRangeHL(hlOptions, spM);
  range_t top_r = validTupleOrthoRangeLH(lhOptions, spM);

  /*
   * Calculate the value of cot(θ) for this middle spacepoint.
   */
  float cotTheta =
      std::max(std::abs(spM.zr()[0] / spM.zr()[1]), lhOptions.cotThetaMax);

  /*
   * Calculate the maximum Δr, given that we have already constrained our
   * search space.
   */
  float deltaRMaxTop = top_r[DimR].max() - spM.zr()[1];
  float deltaRMaxBottom = spM.zr()[1] - bottom_r[DimR].min();

  /*
   * Create the search range for the bottom spacepoint assuming a
   * monotonically increasing z track, by calculating the minimum z value from
   * the cot(θ), and by setting the maximum to the z position of the middle
   * spacepoint - if the z position is higher than the middle point, then it
   * would be a decreasing z track!
   */
  range_t bottom_lh_r = bottom_r;
  bottom_lh_r[DimZ].shrink(spM.zr()[0] - cotTheta * deltaRMaxBottom,
                           spM.zr()[0]);

  /*
   * Calculate the search ranges for the other four sets of points in a
   * similar fashion.
   */
  range_t top_lh_r = top_r;
  top_lh_r[DimZ].shrink(spM.zr()[0], spM.zr()[0] + cotTheta * deltaRMaxTop);

  range_t bottom_hl_r = bottom_r;
  bottom_hl_r[DimZ].shrink(spM.zr()[0],
                           spM.zr()[0] + cotTheta * deltaRMaxBottom);
  range_t top_hl_r = top_r;
  top_hl_r[DimZ].shrink(spM.zr()[0] - cotTheta * deltaRMaxTop, spM.zr()[0]);

  /*
   * Now, we will actually search for the spaces. Remembering that we combine
   * bottom and top candidates for increasing and decreasing tracks
   * separately, we will first check whether both the search ranges for
   * increasing tracks are not degenerate - if they are, we will never find
   * any seeds and we do not need to bother doing the search.
   */
  if (!bottom_lh_r.degenerate() && !top_lh_r.degenerate()) {
    /*
     * Search the trees for points that lie in the given search range.
     */
    m_tree.rangeSearchMapDiscard(
        top_lh_r,
        [&candidates](const Tree::coordinate_t &, const Tree::value_t &top) {
          candidates.top_lh_v.push_back(top);
        });
  }

  /*
   * Perform the same search for candidate bottom spacepoints, but for
   * monotonically decreasing z tracks.
   */
  if (!bottom_hl_r.degenerate() && !top_hl_r.degenerate()) {
    m_tree.rangeSearchMapDiscard(
        top_hl_r,
        [&candidates](const Tree::coordinate_t &, const Tree::value_t &top) {
          candidates.top_hl_v.push_back(top);
        });
  }

  // apply cut on the number of top SP if seedConfirmation is true
  bool search_bot_lh = candidates.top_lh_v.size() >= nTopSeedConf;
  bool search_bot_hl = candidates.top_hl_v.size() >= nTopSeedConf;

  /*
   * Next, we perform a search for bottom candidates in increasing z tracks,
   * which only makes sense if we found any bottom candidates.
   */
  if (!candidates.top_lh_v.empty() && search_bot_lh) {
    m_tree.rangeSearchMapDiscard(
        bottom_lh_r,
        [&candidates](const Tree::coordinate_t &, const Tree::value_t &bottom) {
          candidates.bottom_lh_v.push_back(bottom);
        });
  }

  /*
   * And repeat for the top spacepoints for decreasing z tracks!
   */
  if (!candidates.top_hl_v.empty() && search_bot_hl) {
    m_tree.rangeSearchMapDiscard(
        bottom_hl_r,
        [&candidates](const Tree::coordinate_t &, const Tree::value_t &bottom) {
          candidates.bottom_hl_v.push_back(bottom);
        });
  }
}

CylindricalSpacePointKDTreeBuilder::CylindricalSpacePointKDTreeBuilder(
    std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)) {}

void CylindricalSpacePointKDTreeBuilder::insert(SpacePointIndex index,
                                                float phi, float r, float z) {
  Tree::coordinate_t point{phi, r, z};
  m_points.emplace_back(point, index);
}

void CylindricalSpacePointKDTreeBuilder::extend(
    const SpacePointContainer2::ConstRange &spacePoints) {
  for (const auto &sp : spacePoints) {
    insert(sp);
  }
}

CylindricalSpacePointKDTree CylindricalSpacePointKDTreeBuilder::build() {
  CylindricalSpacePointKDTree result(Tree(std::move(m_points)),
                                     m_logger->clone());
  clear();
  return result;
}

}  // namespace Acts::Experimental
