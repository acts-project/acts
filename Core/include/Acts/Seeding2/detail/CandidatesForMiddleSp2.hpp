// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"

#include <limits>
#include <vector>

namespace Acts::Experimental {

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
  Size nLowQualityCandidates() const { return m_nLow; }

  /// @brief Retrieve the number of High quality candidates
  /// @returns The number of High quality candidates
  Size nHighQualityCandidates() const { return m_nHigh; }

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
  // sizes
  // m_maxSize* is the maximum size of the indices collections. These values
  // are set by the user once
  Size m_maxSizeLow{kNoSize};
  Size m_maxSizeHigh{kNoSize};
  // m_n_* is the current size of the indices collections [0, m_maxSize*).
  // These values are set internally by the class
  Size m_nLow{0};
  Size m_nHigh{0};

  // storage contains the collection of the candidates
  std::vector<TripletCandidate2> m_storage;

  // The following vectors store indexes to elements in the storage
  // They are sorted as a min heap tree, in which
  // Each node is lower than its children
  // Thus, it is guaranteed that the lower elements is at the front
  // Sorting criteria is the seed quality
  //
  // This is in effect faster sorted container - implementation with std::set
  // and std::priority_queue were tried and were found to be slower.

  // list of indexes of candidates with low quality in the storage
  std::vector<Index> m_indicesLow;
  // list of indexes of candidates with high quality in the storage
  std::vector<Index> m_indicesHigh;

  /// @brief Adding a new triplet candidate to the collection, should it satisfy the
  /// selection criteria
  /// @param indices Index collection pointing to the triplet candidates
  /// @param n The current number of stored elements in the container
  /// @param nMax The maximum number of elements that can be stored in the container
  /// @param spB The bottom space point
  /// @param spM The middle space point
  /// @param spT The top space point
  /// @param weight The quality of the triplet candidate
  /// @param zOrigin The z-coordinate of the origin
  /// @param isQuality Whether the triplet candidate is high or low quality
  /// @returns whether the triplet candidate has been added or not to the collection
  bool push(std::vector<Index>& indices, Size& n, Size nMax,
            SpacePointIndex2 spB, SpacePointIndex2 spM, SpacePointIndex2 spT,
            float weight, float zOrigin, bool isQuality);

  /// @brief Check if an element exists in the collection. The element to be checked
  /// is supposed to be in the n position of the collection.
  /// @param n Index of the requested element
  /// @param maxSize Number of elements currently stored in the collection
  /// @returns Whether the element exists
  bool exists(Index n, Size maxSize) const {
    // If the element exists, its index is lower than the current number
    // of stored elements
    return n < maxSize;
  }

  /// @brief Pop an element from a collection. The removal of the element from the collection
  /// does not imply its destruction. In fact, the number of stored elements is
  /// simply diminished by 1. The popped element is technically still available
  /// at the end of the collection.
  /// @param indices Index collection pointing to the triplet candidates
  /// @param currentSize The current number of element stored in the collection. The function will
  /// diminish this value by 1
  void pop(std::vector<Index>& indices, Size& currentSize);

  /// @brief Return the weight for a candidate
  /// @param indices Index collection pointing to the triplet candidates
  /// @param n Index of the element in the collection
  /// @returns The weight of the candidate
  float weight(const std::vector<Index>& indices, Index n) const;

  /// @brief Move an element up in the min heap tree. The function checks whether the element's
  /// weight is lower of its parent's weight. If so, it swaps them. Reiterate
  /// the process until the element is in the correct position on the tree
  /// @param indices Index collection pointing to the triplet candidates
  /// @param n The index of the element to place in the correct position
  void bubbleup(std::vector<Index>& indices, Index n);

  /// @brief Move an element down in the min heap tree. The function checks whether the elements's
  /// weight is lower of its child's weights. If so, it swaps the element with
  /// the child with the lowest weight. Reiterate the process until the element
  /// is in the correct position on the tree
  /// @param indices Index collection pointing to the triplet candidates
  /// @param n The index of the element to place in the correct position
  /// @param actualSize The current number of elements stored in the collection
  void bubbledw(std::vector<Index>& indices, Index n, Size actualSize);

  /// @brief Adding a new triplet candidate to the collection. The function is called after the candidate has satisfied
  /// all the selection criteria
  /// @param indices Index collection pointing to the triplet candidates
  /// @param n Current number of stored elements in the collection
  /// @param nMax The maximum number of elements that can be stored in the collection
  /// @param element The element that must be added to the collection
  void addToCollection(std::vector<Index>& indices, Index& n, Size nMax,
                       const TripletCandidate2& element);
};

}  // namespace Acts::Experimental
