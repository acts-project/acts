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
  
  template<typename container_t>
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
    ContainerType& container();
    const ContainerType& container() const;
    
  private:
    Acts::detail_tc::RefHolder<ContainerType> m_container; 
    IndexType m_index;
  };
  
  // Implementation
  template<typename container_t>
    SpacePointProxy<container_t>::SpacePointProxy(
      typename SpacePointProxy<container_t>::ContainerType& container, 
      typename SpacePointProxy<container_t>::IndexType index)
    : m_container(container), 
    m_index(index)
    {}
  
  template<typename container_t>
    SpacePointProxy<container_t>::SpacePointProxy(
      const typename SpacePointProxy<container_t>::ContainerType& container, 
      typename SpacePointProxy<container_t>::IndexType index)
    : m_container(container), 
    m_index(index)
    {}
  
  template<typename container_t>
    inline typename SpacePointProxy<container_t>::IndexType 
    SpacePointProxy<container_t>::index() const
    { return m_index; }
  
  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::x() const
    { return container().x(m_index); }
  
  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::y() const
    { return container().y(m_index); }
  
  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::z() const
    { return container().z(m_index); }
  
  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::radius() const
    { return container().radius(m_index); }
  
  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::varianceR() const
    { return container().varianceR(m_index); }

  template<typename container_t>
    inline float 
    SpacePointProxy<container_t>::varianceZ() const
    { return container().varianceZ(m_index); }

  template<typename container_t>
    inline typename SpacePointProxy<container_t>::ContainerType& 
    SpacePointProxy<container_t>::container() 
    { return *m_container; }
  
  template<typename container_t>
    inline const typename SpacePointProxy<container_t>::ContainerType& 
    SpacePointProxy<container_t>::container() const
    { return *m_container; }

}
