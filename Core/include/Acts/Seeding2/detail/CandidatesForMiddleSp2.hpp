// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <limits>
#include <vector>

namespace Acts {

/// @brief A description of a triplet candidate.
struct TripletCandidate2 {
  /// @brief Default Constructor
  TripletCandidate2() = default;

  /// @brief constructor
  /// @param b The bottom space point
  /// @param m The middle space point
  /// @param t The top space point
  /// @param w The quality of the candidate
  /// @param z The z coordinate of the origin
  /// @param q Whether the candidate is high or low quality
  TripletCandidate2(SpacePointIndex2 b, SpacePointIndex2 m, SpacePointIndex2 t,
                    float w, float z, bool q)
      : bottom(b), middle(m), top(t), weight(w), zOrigin(z), isQuality(q) {}

  SpacePointIndex2 bottom{};
  SpacePointIndex2 middle{};
  SpacePointIndex2 top{};
  float weight{};
  float zOrigin{};
  bool isQuality{};
};

class CandidatesForMiddleSp2 {
 public:
  using Index = std::uint32_t;
  using Size = std::uint32_t;

  static constexpr Size kNoSize = std::numeric_limits<Size>::max();

  CandidatesForMiddleSp2();

  /// @brief Setting maximum number of candidates to keep
  /// @param nLow Maximum number of candidates in the low-quality collection
  /// @param nHigh Maximum number of candidates in the high-quality collection
  CandidatesForMiddleSp2(Size nLow, Size nHigh);

  Size size() const { return m_storage.size(); }

  /// @brief Clear the internal storage
  void clear();

  /// @brief Retrieve the number of Low quality candidates
  /// @returns The number of Low quality candidates
  Size nLowQualityCandidates() const {
    return static_cast<Size>(m_indicesLow.size());
  }

  /// @brief Retrieve the number of High quality candidates
  /// @returns The number of High quality candidates
  Size nHighQualityCandidates() const {
    return static_cast<Size>(m_indicesHigh.size());
  }

  /// @brief Adding a new triplet candidate to the collection, should it satisfy the
  /// selection criteria
  /// @param spB Bottom space point
  /// @param spM Medium space point
  /// @param spT Top space point
  /// @param weight The quality of the triplet candidate
  /// @param zOrigin The z-coordinate of the origin
  /// @param isQuality Whether the triplet candidate is high or low quality
  /// @returns whether the triplet candidate has been added or not to the collection
  bool push(SpacePointIndex2 spB, SpacePointIndex2 spM, SpacePointIndex2 spT,
            float weight, float zOrigin, bool isQuality);

  /// @brief Retrieve the triplet candidates, the resulting vector is already sorted,
  /// elements with higher quality first
  void toSortedCandidates(std::vector<TripletCandidate2>& output);

 private:
  using WeightIndex = std::pair<float, Index>;
  using Container = std::vector<WeightIndex>;

  static constexpr bool comparator(const WeightIndex& a, const WeightIndex& b) {
    return a.first > b.first;
  }

  // sizes
  // m_maxSize* is the maximum size of the indices collections. These values
  // are set by the user once
  Size m_maxSizeLow{kNoSize};
  Size m_maxSizeHigh{kNoSize};

  // storage contains the collection of the candidates
  std::vector<TripletCandidate2> m_storage;

  Container m_indicesLow;
  Container m_indicesHigh;

  bool push(Container& container, Size nMax, SpacePointIndex2 spB,
            SpacePointIndex2 spM, SpacePointIndex2 spT, float weight,
            float zOrigin, bool isQuality);
};

}  // namespace Acts
