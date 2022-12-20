// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>
#include <tuple>

namespace Acts {

  template<typename external_space_point_t>
  class CandidatesForSpM {

    // Internal collection, we store only the components that are different and useful
    enum InternalComponents : int {
      I_BSP=0,
      I_TSP,
      I_WEIGHT,
      I_ZORIGIN
    };

  public:
    // more complete collection of variables, used by external seeding code
    enum Components : int {
      BSP=0,
      MSP,
      TSP,
      WEIGHT,
      ZORIGIN,
      QUALITY
    };

    using sp_type = external_space_point_t*;
    using value_type = std::tuple<sp_type, sp_type, float, float>;
    using output_type = std::tuple<sp_type, sp_type, sp_type, float, float, bool>;
    static constexpr sp_type default_value = nullptr;
    
    CandidatesForSpM();
    ~CandidatesForSpM() = default;

    void setMaxElements(std::size_t n_low,
			std::size_t n_high);
    void setMediumSp(sp_type);
    void setBottomSp(sp_type);
    const std::vector<value_type>& storage(bool isQuality) const;
    std::vector< output_type > extendedStorage() const;
    const sp_type& spM() const;
    
    void push(sp_type& SpT, float weight, float zOrigin, bool isQuality);
    void clear();
    
  private:
    bool exists(std::size_t, std::size_t) const;

    void pop(std::vector< value_type >&, std::size_t&);
    float top(const std::vector< value_type >&) const;
    float weight(const std::vector< value_type >&, std::size_t) const;

    void bubbleup(std::vector< value_type >&, std::size_t);
    void bubbledw(std::vector< value_type >&, std::size_t, std::size_t);
    
    void addToCollection(std::vector< value_type >&,
			 sp_type& SpB, sp_type& SpT, float weight, float zOrigin,
			 bool isQuality);
    void insertToCollection(std::vector< value_type >&,
			    sp_type& SpB, sp_type& SpT, float weight, float zOrigin,
			    bool isQuality);

  public:
    // sizes
    std::size_t m_max_size_high;
    std::size_t m_max_size_low;
    std::size_t m_n_high;
    std::size_t m_n_low;

    // space points
    sp_type m_SpB;
    sp_type m_SpM;

    // storage
    // These vectors are sorted as a min heap tree
    // Each node is lower then its childs
    // Thus, it is guaranteed that the lower elements is at the front
    // Sorting criteria is the seed quality 

    // storage for candidates with high quality
    std::vector< value_type > m_storage_high;
    // storage for candidates with low quality
    std::vector< value_type > m_storage_low;
  };

  template<typename external_space_point_t>
  inline
  const std::vector<typename CandidatesForSpM<external_space_point_t>::value_type>&
  CandidatesForSpM<external_space_point_t>::storage(bool isQuality) const
  { return isQuality ? m_storage_high : m_storage_low; }

  template<typename external_space_point_t>
  inline
  const typename CandidatesForSpM<external_space_point_t>::sp_type&
  CandidatesForSpM<external_space_point_t>::spM() const
  { return m_SpM; }
  
  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setMaxElements(std::size_t n_low,
								       std::size_t n_high)
  {
    if (m_storage_high.capacity() < n_high) m_storage_high.reserve(n_high);
    if (m_storage_low.capacity() < n_low) m_storage_low.reserve(n_low);
    m_max_size_high = n_high;
    m_max_size_low = n_low;
  }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setMediumSp(typename CandidatesForSpM<external_space_point_t>::sp_type idx)
  { m_SpM = idx; }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::setBottomSp(typename CandidatesForSpM<external_space_point_t>::sp_type idx)
  { m_SpB = idx; }

  template<typename external_space_point_t>
  inline float CandidatesForSpM<external_space_point_t>::top(const std::vector<value_type>& storage) const
  { return weight(storage, 0); }

  template<typename external_space_point_t>
  inline bool CandidatesForSpM<external_space_point_t>::exists(std::size_t n, std::size_t max_size) const
  { return n < max_size; }

  template<typename external_space_point_t>
  inline float CandidatesForSpM<external_space_point_t>::weight(const std::vector<value_type>& storage, std::size_t n) const
  { return std::get<InternalComponents::I_WEIGHT>(storage[n]); }

  template<typename external_space_point_t>
  inline void CandidatesForSpM<external_space_point_t>::clear()
  {
    // do not clear max size, this is set only once
    m_n_high = 0;
    m_n_low = 0;
    // clean fixed space points
    m_SpB = default_value;
    m_SpM = default_value;
    // clean storage
    m_storage_high.clear();
    m_storage_low.clear();
  }
  
}  // namespace Acts

#include "Acts/Seeding/CandidatesForSpM.ipp"
