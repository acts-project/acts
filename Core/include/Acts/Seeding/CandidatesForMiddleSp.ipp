// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::setMaxElements(
    std::size_t n_low, std::size_t n_high) {
  m_max_size_high = n_high;
  m_max_size_low = n_low;

  // protection against default numbers
  // it may cause std::bad_alloc if we don't protect
  if (n_high == std::numeric_limits<std::size_t>::max() or
      n_low == std::numeric_limits<std::size_t>::max()) {
    return;
  }

  // Reserve enough memory for all collections
  m_storage.reserve(n_low + n_high);
  m_indices_high.reserve(n_high);
  m_indices_low.reserve(n_low);
}

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::pop(
    std::vector<std::size_t>& indices, std::size_t& current_size) {
  // Remove the candidate with the lowest weight in the collection
  // By construction, this element is always the first element in the
  // collection.
  // The removal works this way: the first element of the collection
  // is overwritten with the last element of the collection.
  // The number of stored elements (current_size) is lowered by 1
  // The new first element is moved down in the tree to its correct position
  std::swap(indices[0], indices[current_size - 1]);
  bubbledw(indices, 0, --current_size);
}

template <typename external_space_point_t>
inline bool CandidatesForMiddleSp<external_space_point_t>::exists(
    const std::size_t n, const std::size_t max_size) const {
  // If the element exists, its index is lower than the current number
  // of stored elements
  return n < max_size;
}

template <typename external_space_point_t>
inline float CandidatesForMiddleSp<external_space_point_t>::weight(
    const std::vector<std::size_t>& indices, std::size_t n) const {
  // Get the weight of the n-th element
  return m_storage[indices[n]].weight;
}

template <typename external_space_point_t>
inline void CandidatesForMiddleSp<external_space_point_t>::clear() {
  // do not clear max size, this is set only once
  // reset to 0 the number of stored elements
  m_n_high = 0;
  m_n_low = 0;
  // Reset all values to default
  m_storage.clear();
  m_indices_high.clear();
  m_indices_low.clear();
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::push(
    external_space_point_t& SpB, external_space_point_t& SpM,
    external_space_point_t& SpT, float weight, float zOrigin, bool isQuality) {
  // Decide in which collection this candidate may be added to according to the
  // isQuality boolean
  if (isQuality) {
    return push(m_indices_high, m_n_high, m_max_size_high, SpB, SpM, SpT,
                weight, zOrigin, isQuality);
  }
  return push(m_indices_low, m_n_low, m_max_size_low, SpB, SpM, SpT, weight,
              zOrigin, isQuality);
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::push(
    std::vector<std::size_t>& indices, std::size_t& n, const std::size_t n_max,
    external_space_point_t& SpB, external_space_point_t& SpM,
    external_space_point_t& SpT, float weight, float zOrigin, bool isQuality) {
  // If we do not want to store candidates, returns
  if (n_max == 0) {
    return false;
  }

  // if there is still space, add anything
  if (n < n_max) {
    addToCollection(indices, n, n_max,
                    value_type(SpB, SpM, SpT, weight, zOrigin, isQuality));
    return true;
  }

  // if no space, replace one if quality is enough
  // compare to element with lowest weight
  const auto& lowest_weight = this->weight(indices, 0);
  if (weight <= lowest_weight) {
    return false;
  }

  // remove element with lower weight and add this one
  pop(indices, n);
  addToCollection(indices, n, n_max,
                  value_type(SpB, SpM, SpT, weight, zOrigin, isQuality));
  return true;
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::addToCollection(
    std::vector<std::size_t>& indices, std::size_t& n, const std::size_t n_max,
    value_type&& element) {
  // adds elements to the end of the collection
  if (indices.size() == n_max) {
    m_storage[indices[n]] = std::move(element);
  } else {
    m_storage.push_back(std::move(element));
    indices.push_back(m_storage.size() - 1);
  }
  // Move the added element up in the tree to its correct position
  bubbleup(indices, n++);
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::bubbledw(
    std::vector<std::size_t>& indices, std::size_t n, std::size_t actual_size) {
  while (n < actual_size) {
    // The collection of indexes are sorted as min heap trees
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    float current = weight(indices, n);
    std::size_t left_child = 2 * n + 1;
    std::size_t right_child = 2 * n + 2;

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
    if (not exists(left_child, actual_size)) {
      break;
    }

    // At least one of the child is present. Left child for sure, right child we
    // have to check. We take the lowest weight of the children. By default this
    // is the weight of the left child, and we then check for the right child

    float weight_left_child = weight(indices, left_child);

    std::size_t selected_child = left_child;
    float selected_weight = weight_left_child;

    // Check which child has the lower weight
    if (exists(right_child, actual_size)) {
      float weight_right_child = weight(indices, right_child);
      if (weight_right_child <= weight_left_child) {
        selected_child = right_child;
        selected_weight = weight_right_child;
      }
    }

    // At this point we have the minimum weight of the children
    // We can compare this to the current weight
    // If weight of the children is higher we stop
    if (selected_weight >= current) {
      break;
    }

    // swap and repeat the process
    std::swap(indices[n], indices[selected_child]);
    n = selected_child;
  }  // while loop
}

template <typename external_space_point_t>
void CandidatesForMiddleSp<external_space_point_t>::bubbleup(
    std::vector<std::size_t>& indices, std::size_t n) {
  while (n != 0) {
    // The collection of indexes are sorted as min heap trees
    // parent: (n - 1) / 2;
    // this works because it is an integer operation
    std::size_t parent_idx = (n - 1) / 2;

    float weight_current = weight(indices, n);
    float weight_parent = weight(indices, parent_idx);

    // If weight of the parent is lower than this one, we stop
    if (weight_parent <= weight_current) {
      break;
    }

    // swap and repeat the process
    std::swap(indices[n], indices[parent_idx]);
    n = parent_idx;
  }
}

template <typename external_space_point_t>
std::vector<typename CandidatesForMiddleSp<external_space_point_t>::value_type>
CandidatesForMiddleSp<external_space_point_t>::storage() {
  // this will retrieve the entire storage
  // the resulting vector is already sorted from high to low quality
  std::vector<value_type> output(m_n_high + m_n_low);
  std::size_t out_idx = output.size() - 1;

  // rely on the fact that m_indices_* are both min heap trees
  // Sorting comes naturally by popping elements one by one and
  // placing this element at the end of the output vector
  while (m_n_high != 0 or m_n_low != 0) {
    // no entries in collection high, we attach the entire low collection
    if (m_n_high == 0) {
      std::size_t idx = m_n_low;
      for (std::size_t i(0); i < idx; i++) {
        output[out_idx--] = std::move(m_storage[m_indices_low[0]]);
        pop(m_indices_low, m_n_low);
      }
      break;
    }

    // no entries in collection low, we attach the entire high collection
    if (m_n_low == 0) {
      std::size_t idx = m_n_high;
      for (std::size_t i(0); i < idx; i++) {
        output[out_idx--] = std::move(m_storage[m_indices_high[0]]);
        pop(m_indices_high, m_n_high);
      }
      break;
    }

    // Both have entries, get the minimum
    if (descendingByQuality(m_storage[m_indices_low[0]],
                            m_storage[m_indices_high[0]])) {
      output[out_idx--] = std::move(m_storage[m_indices_high[0]]);
      pop(m_indices_high, m_n_high);
    } else {
      output[out_idx--] = std::move(m_storage[m_indices_low[0]]);
      pop(m_indices_low, m_n_low);
    }

  }  // while loop

  clear();
  return output;
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::descendingByQuality(
    const value_type& i1, const value_type& i2) {
  if (i1.weight != i2.weight) {
    return i1.weight > i2.weight;
  }

  // This is for the case when the weights from different seeds
  // are same. This makes cpu & cuda results same

  const auto& bottom_l1 = i1.bottom;
  const auto& middle_l1 = i1.middle;
  const auto& top_l1 = i1.top;

  const auto& bottom_l2 = i2.bottom;
  const auto& middle_l2 = i2.middle;
  const auto& top_l2 = i2.top;

  float seed1_sum = 0.;
  float seed2_sum = 0.;

  seed1_sum +=
      bottom_l1->y() * bottom_l1->y() + bottom_l1->z() * bottom_l1->z();
  seed1_sum +=
      middle_l1->y() * middle_l1->y() + middle_l1->z() * middle_l1->z();
  seed1_sum += top_l1->y() * top_l1->y() + top_l1->z() * top_l1->z();

  seed2_sum +=
      bottom_l2->y() * bottom_l2->y() + bottom_l2->z() * bottom_l2->z();
  seed2_sum +=
      middle_l2->y() * middle_l2->y() + middle_l2->z() * middle_l2->z();
  seed2_sum += top_l2->y() * top_l2->y() + top_l2->z() * top_l2->z();

  return seed1_sum > seed2_sum;
}

template <typename external_space_point_t>
bool CandidatesForMiddleSp<external_space_point_t>::ascendingByQuality(
    const value_type& i1, const value_type& i2) {
  return not descendingByQuality(i1, i2);
}

}  // namespace Acts
