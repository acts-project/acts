// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

  template<typename external_space_point_t>
  CandidatesForSpM<external_space_point_t>::CandidatesForSpM()
    : m_max_size_high(0),
      m_max_size_low(0),
      m_n_high(0),
      m_n_low(0),	
      m_SpB(CandidatesForSpM<external_space_point_t>::default_value),
      m_SpM(CandidatesForSpM<external_space_point_t>::default_value)
  {}

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::push(typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
			      float weight, float zOrigin,
			      bool isQuality)
  {
    auto& storage = isQuality ? m_storage_high : m_storage_low;
    const std::size_t& current_max_size = isQuality ? m_max_size_high : m_max_size_low;
    std::size_t& current_size = isQuality ? m_n_high : m_n_low;

    // if there is still space, add anything
    if (current_size < current_max_size) {
      addToCollection(storage,
		      m_SpB, SpT, weight, zOrigin,
		      isQuality);
      return;
    }

    // if no space, replace one if quality is enough
    // compare to element with lower weight
    const auto& lower_weight = top(storage);
    if (weight <= lower_weight)
      return;
    
    // remove element with lower weight and add this one  
    pop(storage, current_size);
    insertToCollection(storage,
    		       m_SpB, SpT, weight, zOrigin,
		       isQuality);
  }

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::addToCollection(std::vector< value_type >& storage,
       					 typename CandidatesForSpM<external_space_point_t>::sp_type& SpB,
					 typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
					 float weight, float zOrigin,
					 bool isQuality)
  {
    // adds elements to the end of the collection
    // function called when space in storage is not full
    auto toAdd = std::make_tuple(SpB, SpT, weight, zOrigin);
    storage.push_back( toAdd );
    std::size_t& added_index = isQuality ? m_n_high : m_n_low;
    bubbleup(storage, added_index);
    ++added_index;
}  

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::insertToCollection(std::vector< value_type >& storage,
                                            typename CandidatesForSpM<external_space_point_t>::sp_type& SpB,
					    typename CandidatesForSpM<external_space_point_t>::sp_type& SpT,
					    float weight, float zOrigin,
					    bool isQuality)
  {
    // inserts elements to the end of the collection
    // function called when space in storage is full
    // before this a pop is called
    auto toAdd = std::make_tuple(SpB, SpT, weight, zOrigin);
    std::size_t& added_index = isQuality ? m_n_high : m_n_low;
    storage[added_index] = toAdd;
    bubbleup(storage, added_index);
    ++added_index;
}

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::bubbledw(std::vector< value_type >& storage,
       							  std::size_t n,
							  std::size_t actual_size)
  {
    // left child : 2 * n + 1
    // right child: 2 * n + 2
    float current = weight(storage, n);
    std::size_t left_child = 2 * n + 1;
    std::size_t right_child = 2 * n + 2;
    
    // no left child, we stop
    if (not exists(left_child, actual_size)) return;

    float weight_left_child = weight(storage, left_child);
    
    // no right child, left wins
    if (not exists(right_child, actual_size)) {
      if (weight_left_child < current) {
	std::swap(storage[n], storage[left_child]);
	return bubbledw(storage, left_child, actual_size);
      }
    }

    float weight_right_child = weight(storage, right_child);

    // both childs
    // left is smaller
    if (weight_left_child < weight_right_child) {
      if (weight_left_child < current) {
	std::swap(storage[n], storage[left_child]);
	return bubbledw(storage, left_child, actual_size);
      }
    }
    // right is smaller
    if (weight_right_child < current) {
      std::swap(storage[n], storage[right_child]);
      return bubbledw(storage, right_child, actual_size);
    }
    
  }

  template<typename external_space_point_t>	
  void CandidatesForSpM<external_space_point_t>::bubbleup(std::vector< value_type >& storage,
       							  std::size_t n)
  {
    if (n == 0) return;
    
    // parent: (n - 1) / 2;
    // this works because it is an integer operation
    std::size_t parent_idx = (n - 1) / 2;

    float weight_current = weight(storage, n);
    float weight_parent = weight(storage, parent_idx);

    if (weight_parent <= weight_current)
      { return; }
    
    std::swap(storage[n], storage[parent_idx]);
    bubbleup(storage, parent_idx);
  }

  template<typename external_space_point_t>
  void CandidatesForSpM<external_space_point_t>::pop(std::vector< value_type >& storage,
       						     std::size_t& current_size)
  {
    storage[0] = storage[current_size - 1];
    --current_size;
    bubbledw(storage, 0, current_size);
  }

  template<typename external_space_point_t>
  std::vector< typename CandidatesForSpM<external_space_point_t>::output_type >
  CandidatesForSpM<external_space_point_t>::extendedStorage() const
  {
    // this will retrieve the entire storage, first high and then low quality
    // the resulting vector is not sorted!
    std::vector< output_type > output;
    output.reserve(m_n_high + m_n_low);

    for (std::size_t idx(0); idx < m_n_high; idx++) {
    	const auto& [bottom, top, weight, zOrigin] = m_storage_high[idx];
	output.emplace_back( bottom, m_SpM, top, weight, zOrigin, true );
    }
	
    for (std::size_t idx(0); idx < m_n_low; idx++) {
       const auto& [bottom, top, weight, zOrigin] = m_storage_low[idx];
       output.emplace_back( bottom, m_SpM, top, weight, zOrigin, false );
    }

    // sort output according to weight and sps
    // should we collect inputs according to this criterion instead?
    std::sort(output.begin(), output.end(),
            [] (const auto& i1, const auto& i2) -> bool
        {
	    const auto& [bottom_l1, medium_l1, top_l1, weight_l1, zOrigin_l1, isQuality_l1] = i1;
	    const auto& [bottom_l2, medium_l2, top_l2, weight_l2, zOrigin_l2, isQuality_l2] = i2;

	    if (weight_l1 != weight_l2)
	       return weight_l1 > weight_l2;

	    // This is for the case when the weights from different seeds
	    // are same. This makes cpu & cuda results same

	    // medium is the same for all candidates
	    float sum_medium = medium_l1->y() * medium_l1->y() + medium_l1->z() * medium_l1->z();

	    float seed1_sum = sum_medium;
	    float seed2_sum = sum_medium;

	    seed1_sum += bottom_l1->y() * bottom_l1->y() + bottom_l1->z() * bottom_l1->z();
	    seed1_sum += top_l1->y() * top_l1->y() + top_l1->z() * top_l1->z();

	    seed2_sum += bottom_l2->y() * bottom_l2->y() + bottom_l2->z() * bottom_l2->z();
	    seed2_sum += top_l2->y() * top_l2->y() + top_l2->z() * top_l2->z();

	    return seed1_sum > seed2_sum;
      });

    return output;
  }
  
} //namespace

