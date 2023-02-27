// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Holders.hpp"
#include <type_traits>

namespace Acts {
  
  template<typename container_t,
	   bool read_only>
  class SpacePointProxy {
  public:
    using ContainerType = container_t;
    using IndexType = typename ContainerType::IndexType;
    
  public:
    // Never take the ownership of the container
    SpacePointProxy(ContainerType&& container, IndexType index) = delete;
    // Only get the reference
    SpacePointProxy(ContainerType& container, IndexType index);
    SpacePointProxy(const ContainerType& container, IndexType index);

    // copy and move operations are defaults
    
    IndexType index() const;
    float x() const;
    float y() const;
    float z() const;
    float radius() const;
    float varianceR() const;
    float varianceZ() const;
    
    // Add component methods for additional quantities
    
  private:
    template<bool RO = read_only, std::enable_if_t<!RO, bool> = true>
    ContainerType& container() { return *m_container; }
    
    const ContainerType& container() const;
    
  private:
    Acts::detail_tc::RefHolder<ContainerType> m_container; 
    IndexType m_index;
  };
  
  // Implementation
  template<typename container_t, bool read_only>
  SpacePointProxy<container_t, read_only>::SpacePointProxy(
      typename SpacePointProxy<container_t, read_only>::ContainerType& container, 
      typename SpacePointProxy<container_t, read_only>::IndexType index)
    : m_container(container), 
    m_index(index)
    {}
  
  template<typename container_t, bool read_only>
  SpacePointProxy<container_t, read_only>::SpacePointProxy(
      const typename SpacePointProxy<container_t, read_only>::ContainerType& container, 
      typename SpacePointProxy<container_t, read_only>::IndexType index)
    : m_container(container), 
    m_index(index)
    {}
  
  template<typename container_t, bool read_only>
  inline typename SpacePointProxy<container_t, read_only>::IndexType 
  SpacePointProxy<container_t, read_only>::index() const
    { return m_index; }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::x() const
  { return container().x(m_index); }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::y() const
  { return container().y(m_index); }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::z() const
  { return container().z(m_index); }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::radius() const
  { return container().radius(m_index); }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::varianceR() const
  { return container().varianceR(m_index); }
  
  template<typename container_t, bool read_only>
  inline float 
  SpacePointProxy<container_t, read_only>::varianceZ() const
  { return container().varianceZ(m_index); }
  
  template<typename container_t, bool read_only>
  inline const typename SpacePointProxy<container_t, read_only>::ContainerType& 
  SpacePointProxy<container_t, read_only>::container() const
  { return *m_container; }

}
