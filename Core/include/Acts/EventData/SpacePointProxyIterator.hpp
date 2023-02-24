// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Holders.hpp"
#include "Acts/EventData/SpacePointProxy.hpp"

namespace Acts {
  
  template<typename container_t,
	   //	   typename proxy_t, // can it be in container specifications?
	   bool read_only = true>
  class SpacePointProxyIterator {
  public:
    using ContainerType = typename std::conditional<read_only, const container_t, container_t>::type;
    using ProxyType = typename container_t::SpacePointProxyType; 
    using ConstProxyType = const ProxyType;
    using IndexType = typename container_t::IndexType;
    
    using iterator_category = std::random_access_iterator_tag;
    using value_type = typename std::conditional<read_only, ConstProxyType, ProxyType>::type; 
    using difference_type = std::ptrdiff_t;
    using pointer = value_type*;
    using reference = value_type&;
    
    // Constructors  
    SpacePointProxyIterator(ContainerType& container, IndexType index);   
    
    SpacePointProxyIterator& operator++();
    SpacePointProxyIterator& operator--();  
    SpacePointProxyIterator operator++(int);
    SpacePointProxyIterator operator--(int);
    
    bool operator==(const SpacePointProxyIterator& other) const;
    bool operator!=(const SpacePointProxyIterator& other) const;
    bool operator<(const SpacePointProxyIterator& other) const;
    bool operator>(const SpacePointProxyIterator& other) const;
    bool operator<=(const SpacePointProxyIterator& other) const;
    bool operator>=(const SpacePointProxyIterator& other) const;
    
    SpacePointProxyIterator& operator+=(IndexType offset);
    SpacePointProxyIterator& operator-=(IndexType offset);
    
    SpacePointProxyIterator operator+(IndexType offset) const;
    SpacePointProxyIterator operator-(IndexType offset) const;
    
    // returning pointer here since seeding always assumes a vector of pointers as input...
    // not sure about this tbh 
    ConstProxyType operator*() const;
    
    template<bool RO = read_only, typename = std::enable_if_t<!RO>>
    ProxyType operator*();
    
  private:
    Acts::detail_tc::RefHolder<ContainerType> m_container;
    IndexType m_index;
  };
  
  // Implementation
  template<typename container_t, bool read_only>
    SpacePointProxyIterator<container_t, read_only>::SpacePointProxyIterator(
      typename SpacePointProxyIterator<container_t, read_only>::ContainerType& container,
      typename SpacePointProxyIterator<container_t, read_only>::IndexType index)
    : m_container(container), 
    m_index(index)
    {}
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>& 
    SpacePointProxyIterator<container_t, read_only>::operator++() 
    {
      ++m_index;
      return *this;
    }
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>&
    SpacePointProxyIterator<container_t, read_only>::operator--() 
    {
      --m_index;
      return *this;
    }
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>
    SpacePointProxyIterator<container_t, read_only>::operator++(int) 
    {
      SpacePointProxyIterator other(*this);
      ++m_index;
      return other;
    }
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>
    SpacePointProxyIterator<container_t, read_only>::operator--(int) 
    {
      SpacePointProxyIterator other(*this);
      --m_index;
      return other;
    }
  
  template<typename container_t, bool read_only>
    inline bool
    SpacePointProxyIterator<container_t, read_only>::operator==(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return m_container.ptr == other.m_container.ptr and
      m_index == other.m_index;;
    }
  
  template<typename container_t, bool read_only>
    inline bool 
    SpacePointProxyIterator<container_t, read_only>::operator!=(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return not (*this == other) ;
    }

  template<typename container_t, bool read_only>  
    inline bool 
    SpacePointProxyIterator<container_t, read_only>::operator<(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return m_index < other.m_index;
    }

  template<typename container_t, bool read_only>
    inline bool 
    SpacePointProxyIterator<container_t, read_only>::operator>(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return m_index > other.m_index;
    }

  template<typename container_t, bool read_only>  
    inline bool 
    SpacePointProxyIterator<container_t, read_only>::operator<=(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return m_index <= other.m_index;
    }

  template<typename container_t, bool read_only>  
    inline bool 
    SpacePointProxyIterator<container_t, read_only>::operator>=(
      const SpacePointProxyIterator<container_t, read_only>& other) const 
    {
      return m_index >= other.m_index;
    }

  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>& 
    SpacePointProxyIterator<container_t, read_only>::operator+=(
      typename SpacePointProxyIterator<container_t, read_only>::IndexType offset) 
    {
      m_index += offset;
      return *this;
    }
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>& 
    SpacePointProxyIterator<container_t, read_only>::operator-=(
      typename SpacePointProxyIterator<container_t, read_only>::IndexType offset) 
    {
      m_index -= offset;
      return *this;
    }

  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>
    SpacePointProxyIterator<container_t, read_only>::operator+(
      typename SpacePointProxyIterator<container_t, read_only>::IndexType offset) const 
    {
      return SpacePointProxyIterator(*m_container, m_index + offset);
    }
  
  template<typename container_t, bool read_only>
    inline SpacePointProxyIterator<container_t, read_only>
    SpacePointProxyIterator<container_t, read_only>::operator-(
      typename SpacePointProxyIterator<container_t, read_only>::IndexType offset) const 
    {
      return SpacePointProxyIterator(*m_container, m_index - offset);
    }

  
  template<typename container_t, bool read_only>
    template<bool, typename>
    inline typename SpacePointProxyIterator<container_t, read_only>::ProxyType
    SpacePointProxyIterator<container_t, read_only>::operator*() 
    { return m_container->get(m_index); }
  
  template<typename container_t, bool read_only>
    inline typename SpacePointProxyIterator<container_t, read_only>::ConstProxyType
    SpacePointProxyIterator<container_t, read_only>::operator*() const
    { return m_container->get(m_index); }
  
}
