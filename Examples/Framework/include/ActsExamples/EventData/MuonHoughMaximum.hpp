// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/EventData/MuonSpacePoint.hpp"

namespace ActsExamples {
/** @brief EDM class to store the result from the HoughTransform pattern finding step. The class takes the
 *         local line parameters of the transform in the precision direction
 * (tanAlpha, interceptY) if the transform, is performed to find an parameter
 * estimate in the bending direction or (tanAlpha, interceptY, tanBeta,
 * interceptX) if a second hough transform is performed on top of the first
 * transform to find the parameters in the complementary direction. Further, the
 * MuonHoughMaximum stores the list of all hits associated with the maximum. */
class MuonHoughMaximum {
 public:
  using HitVec = std::vector<const MuonSpacePoint*>;
  /** @brief Constructor taking the estimated hough parameters and the associated hits
   *  @param tanBeta: Slope of the estimated line in precision direction
   *  @param interceptY: Intercept of the line along the precision direction
   *  @param assocHits: List of hits associated with these parameters */
  MuonHoughMaximum(const double tanBeta, const double interceptY,
                   const HitVec& assocHits)
      : m_tanBeta{tanBeta}, m_interceptY{interceptY}, m_hits{assocHits} {}
  /** @brief Constructor taking the estimated hough parameters from a complementary hough transform
   *        & the associated hits.
   *  @param tanAlpha: Slope of the estimate line in the non-bending direction
   *  @param interceptX: Intercept of the line along the non-bending direction
   *  @param tanBeta: Slope of the estimated line in precision direction
   *  @param interceptY: Intercept of the line along the precision direction
   *  @param assocHits: List of hits associated with these parameters */
  MuonHoughMaximum(const double tanAlpha, const double interceptX,
                   const double tanBeta, const double interceptY,
                   const HitVec& assocHits)
      : m_tanAlpha{tanAlpha},
        m_interceptX{interceptX},
        m_tanBeta{tanBeta},
        m_interceptY{interceptY},
        m_hits{assocHits} {}
  /** @brief Return the slope along the non-bending direction */
  double tanAlpha() const { return m_tanAlpha; }
  /** @brief Return slope along the bending direction */
  double tanBeta() const { return m_tanBeta; }
  /** @brief Return the intercept in the non-precision plane */
  double interceptX() const { return m_interceptX; }
  /** @brief Return the intercept in the precision plane */
  double interceptY() const { return m_interceptY; }
  /** @brief Return the associated hits */
  const HitVec& hits() const { return m_hits; }

 private:
  double m_tanAlpha{0.};
  double m_interceptX{0.};
  double m_tanBeta{0.};
  double m_interceptY{0.};
  HitVec m_hits{};
};

using MuonHoughMaxContainer = std::vector<MuonHoughMaximum>;
}  // namespace ActsExamples
