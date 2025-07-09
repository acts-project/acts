// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/detail/CandidatesForMiddleSp2.hpp"

#include <limits>

namespace Acts::Experimental {

void CandidatesForMiddleSp2::setMaxElements(Size nLow, Size nHigh) {
  m_maxSizeHigh = nHigh;
  m_maxSizeLow = nLow;

  // protection against default numbers
  // it may cause std::bad_alloc if we don't protect
  if (nHigh == std::numeric_limits<Size>::max() ||
      nLow == std::numeric_limits<Size>::max()) {
    return;
  }

  // Reserve enough memory for all collections
  m_storage.reserve(nLow + nHigh);
  m_indicesHigh.reserve(nHigh);
  m_indicesLow.reserve(nLow);
}

void CandidatesForMiddleSp2::pop(std::vector<Index>& indices,
                                 Size& currentSize) {
  // Remove the candidate with the lowest weight in the collection
  // By construction, this element is always the first element in the
  // collection.
  // The removal works this way: the first element of the collection
  // is overwritten with the last element of the collection.
  // The number of stored elements (currentSize) is lowered by 1
  // The new first element is moved down in the tree to its correct position
  std::swap(indices[0], indices[currentSize - 1]);
  bubbledw(indices, 0, --currentSize);
}

float CandidatesForMiddleSp2::weight(const std::vector<Index>& indices,
                                     const Index n) const {
  // Get the weight of the n-th element
  return m_storage[indices[n]].weight;
}

void CandidatesForMiddleSp2::clear() {
  // do not clear max size, this is set only once
  // reset to 0 the number of stored elements
  m_nHigh = 0;
  m_nLow = 0;
  // Reset all values to default
  m_storage.clear();
  m_indicesHigh.clear();
  m_indicesLow.clear();
}

bool CandidatesForMiddleSp2::push(SpacePointIndex2 spB, SpacePointIndex2 spM,
                                  SpacePointIndex2 spT, float weight,
                                  float zOrigin, bool isQuality) {
  // Decide in which collection this candidate may be added to according to the
  // isQuality boolean
  if (isQuality) {
    return push(m_indicesHigh, m_nHigh, m_maxSizeHigh, spB, spM, spT, weight,
                zOrigin, isQuality);
  }
  return push(m_indicesLow, m_nLow, m_maxSizeLow, spB, spM, spT, weight,
              zOrigin, isQuality);
}

bool CandidatesForMiddleSp2::push(std::vector<Index>& indices, Index& n,
                                  const Size nMax, const Index spB,
                                  const Index spM, const Index spT,
                                  const float weight, const float zOrigin,
                                  const bool isQuality) {
  // If we do not want to store candidates, returns
  if (nMax == 0) {
    return false;
  }

  // if there is still space, add anything
  if (n < nMax) {
    addToCollection(
        indices, n, nMax,
        TripletCandidate2(spB, spM, spT, weight, zOrigin, isQuality));
    return true;
  }

  // if no space, replace one if quality is enough
  // compare to element with lowest weight
  if (float lowestWeight = this->weight(indices, 0); weight <= lowestWeight) {
    return false;
  }

  // remove element with lower weight and add this one
  pop(indices, n);
  addToCollection(indices, n, nMax,
                  TripletCandidate2(spB, spM, spT, weight, zOrigin, isQuality));
  return true;
}

void CandidatesForMiddleSp2::addToCollection(std::vector<Index>& indices,
                                             Index& n, const Size nMax,
                                             const TripletCandidate2& element) {
  // adds elements to the end of the collection
  if (indices.size() == nMax) {
    m_storage[indices[n]] = element;
  } else {
    m_storage.push_back(element);
    indices.push_back(m_storage.size() - 1);
  }
  // Move the added element up in the tree to its correct position
  bubbleup(indices, n++);
}

void CandidatesForMiddleSp2::bubbledw(std::vector<Index>& indices, Index n,
                                      const Size actualSize) {
  while (n < actualSize) {
    // The collection of indexes are sorted as min heap trees
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    const float current = weight(indices, n);
    const Index leftChild = 2 * n + 1;
    const Index rightChild = 2 * n + 2;

    // We have to move the current node down the tree to its correct position.
    // This is done by comparing its weight with the weights of its two
    // children. Few things can happen:
    //   - there are no children
    //   - the current weight is lower than the weight of the children
    //   - at least one of the children has a lower weight
    // In the first two cases we stop, since we are already in the correct
    // position

    // if there is no left child, that also means no right child is present.
    // We do nothing
    if (!exists(leftChild, actualSize)) {
      break;
    }

    // At least one of the children is present. Left child for sure, right child
    // we have to check. We take the lowest weight of the children. By default
    // this is the weight of the left child, and we then check for the right
    // child

    const float weightLeftChild = weight(indices, leftChild);

    Index selectedChild = leftChild;
    float selectedWeight = weightLeftChild;

    // Check which child has the lower weight
    if (exists(rightChild, actualSize)) {
      float weightRightChild = weight(indices, rightChild);
      if (weightRightChild <= weightLeftChild) {
        selectedChild = rightChild;
        selectedWeight = weightRightChild;
      }
    }

    // At this point we have the minimum weight of the children
    // We can compare this to the current weight
    // If weight of the children is higher we stop
    if (selectedWeight >= current) {
      break;
    }

    // swap and repeat the process
    std::swap(indices[n], indices[selectedChild]);
    n = selectedChild;
  }  // while loop
}

void CandidatesForMiddleSp2::bubbleup(std::vector<Index>& indices, Index n) {
  while (n != 0) {
    // The collection of indexes are sorted as min heap trees
    // parent: (n - 1) / 2;
    // this works because it is an integer operation
    const Index parentIdx = (n - 1) / 2;

    const float weightCurrent = weight(indices, n);
    // If weight of the parent is lower than this one, we stop
    if (float weightParent = weight(indices, parentIdx);
        weightParent <= weightCurrent) {
      break;
    }

    // swap and repeat the process
    std::swap(indices[n], indices[parentIdx]);
    n = parentIdx;
  }
}

void CandidatesForMiddleSp2::toSortedCandidates(
    const SpacePointContainer2& spacePoints,
    std::vector<TripletCandidate2>& output) {
  // this will retrieve the entire storage
  // the resulting vector is already sorted from high to low quality
  output.resize(m_nHigh + m_nLow);
  Index outIdx = output.size() - 1;

  // rely on the fact that m_indices* are both min heap trees
  // Sorting comes naturally by popping elements one by one and
  // placing this element at the end of the output vector
  while (m_nHigh != 0 || m_nLow != 0) {
    // no entries in collection high, we attach the entire low collection
    if (m_nHigh == 0) {
      Index idx = m_nLow;
      for (Index i = 0; i < idx; i++) {
        output[outIdx] = m_storage[m_indicesLow[0]];
        pop(m_indicesLow, m_nLow);
        outIdx--;
      }
      break;
    }

    // no entries in collection low, we attach the entire high collection
    if (m_nLow == 0) {
      Index idx = m_nHigh;
      for (Index i = 0; i < idx; i++) {
        output[outIdx] = m_storage[m_indicesHigh[0]];
        pop(m_indicesHigh, m_nHigh);
        outIdx--;
      }
      break;
    }

    // Both have entries, get the minimum
    if (descendingByQuality(spacePoints, m_storage[m_indicesLow[0]],
                            m_storage[m_indicesHigh[0]])) {
      output[outIdx--] = m_storage[m_indicesHigh[0]];
      pop(m_indicesHigh, m_nHigh);
    } else {
      output[outIdx--] = m_storage[m_indicesLow[0]];
      pop(m_indicesLow, m_nLow);
    }
  }  // while loop

  clear();
}

bool CandidatesForMiddleSp2::descendingByQuality(
    const SpacePointContainer2& spacePoints, const TripletCandidate2& i1,
    const TripletCandidate2& i2) {
  if (i1.weight != i2.weight) {
    return i1.weight > i2.weight;
  }

  // This is for the case when the weights from different seeds
  // are same. This makes cpu & cuda results same

  const auto bottomL1 = spacePoints[i1.bottom];
  const auto middleL1 = spacePoints[i1.middle];
  const auto topL1 = spacePoints[i1.top];

  const auto bottomL2 = spacePoints[i2.bottom];
  const auto middleL2 = spacePoints[i2.middle];
  const auto topL2 = spacePoints[i2.top];

  float seed1_sum = 0.;
  float seed2_sum = 0.;

  seed1_sum += bottomL1.y() * bottomL1.y() + bottomL1.z() * bottomL1.z();
  seed1_sum += middleL1.y() * middleL1.y() + middleL1.z() * middleL1.z();
  seed1_sum += topL1.y() * topL1.y() + topL1.z() * topL1.z();

  seed2_sum += bottomL2.y() * bottomL2.y() + bottomL2.z() * bottomL2.z();
  seed2_sum += middleL2.y() * middleL2.y() + middleL2.z() * middleL2.z();
  seed2_sum += topL2.y() * topL2.y() + topL2.z() * topL2.z();

  return seed1_sum > seed2_sum;
}

bool CandidatesForMiddleSp2::ascendingByQuality(
    const SpacePointContainer2& spacePoints, const TripletCandidate2& i1,
    const TripletCandidate2& i2) {
  return !descendingByQuality(spacePoints, i1, i2);
}

}  // namespace Acts::Experimental
