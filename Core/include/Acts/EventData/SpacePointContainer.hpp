// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointProxy.hpp" 
#include "Acts/EventData/SpacePointProxyIterator.hpp"
#include "Acts/EventData/Utils.hpp"

#include <vector>

namespace Acts {
  // add read and write here?
  template<typename container_t,
	   template <typename> class holder_t>
  class SpacePointContainer {
  public:
    friend class Acts::SpacePointProxy<Acts::SpacePointContainer<container_t, holder_t>>;
    friend class Acts::SpacePointProxyIterator<Acts::SpacePointContainer<container_t, holder_t>, 
					       false>;
    friend class Acts::SpacePointProxyIterator<Acts::SpacePointContainer<container_t, holder_t>,
					       true>;

    static constexpr bool ReadOnly = true;

    using IndexType = typename container_t::IndexType;
    using SpacePointProxyType = Acts::SpacePointProxy<Acts::SpacePointContainer<container_t, holder_t>>;
    using iterator = Acts::SpacePointProxyIterator<Acts::SpacePointContainer<container_t, holder_t>,
						   false>;
    using const_iterator = Acts::SpacePointProxyIterator<Acts::SpacePointContainer<container_t, holder_t>,
							 true>;
    
  public:    
    // Constructors
    // It makes sense to support both options of
    // taking or not the ownership

    // Do not take ownership
    // Activate only if holder_t is RefHolder
    template<template <typename> class H = holder_t,
             std::enable_if_t<
               Acts::detail_tc::is_same_template<H, Acts::detail_tc::RefHolder>::value
               >* = nullptr>
    SpacePointContainer(container_t& container)
      : m_container( container )
    {
      IndexType n = size();
      m_proxies.reserve(n);
      for (IndexType i(0); i<n; i++)
        m_proxies.emplace_back( *this, i );
    }
    
    // Take the ownership
    // Activate only if holder_t is ValueHolder
    template<template <typename> class H = holder_t,
             std::enable_if_t<
               Acts::detail_tc::is_same_template<H, Acts::detail_tc::ValueHolder>::value
               >* = nullptr>
    SpacePointContainer(container_t&& container)
      :	m_container( container )
    {
      IndexType n = size();
      m_proxies.reserve(n);
      for (IndexType i(0); i<n; i++)
        m_proxies.emplace_back( *this, i );
    }
    
    // If we take ownership, forbid copy operations
    // Need to define copy operations only if holder_t is RefHolder !!!
    template<template <typename> class H = holder_t,
	     std::enable_if_t<
	       Acts::detail_tc::is_same_template<H, Acts::detail_tc::RefHolder>::value
	       >* = nullptr>
    SpacePointContainer(SpacePointContainer& other)
      : m_container( *m_container.ptr ),
	m_proxies( other.m_proxies.begin(), other.m_proxies.end() )
    {}

    template<template <typename> class H = holder_t,
	     std::enable_if_t<Acts::detail_tc::is_same_template<H, Acts::detail_tc::RefHolder>::value, bool> = true>
    SpacePointContainer& operator=(SpacePointContainer& other)
    {
      m_container.ptr = other.m_container.ptr;
      m_proxies.insert( m_proxies.end(), other.m_proxies.begin(), other.m_proxies.end() );
      return *this;
    }
    
    // move operations
    SpacePointContainer(SpacePointContainer&& other) noexcept
      : m_container( std::exchange( other.m_container.ptr, nullptr) ),
	m_proxies( std::move(other.m_proxies) )
    {}

    SpacePointContainer& operator=(SpacePointContainer&& other) noexcept
    {
      m_container = std::exchange( other.m_container.ptr, nullptr);
      m_proxies = std::move( other.m_proxies );
      return *this;
    }
    
    // Destructor
    ~SpacePointContainer() = default;
    
    IndexType size() const;

    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;

    SpacePointProxyType& get(IndexType n);
    const SpacePointProxyType& get(IndexType n) const;

    // do these need to be private or public?
    const container_t& container() const;

  private:
    float x(IndexType n) const;
    float y(IndexType n) const;
    float z(IndexType n) const;
    float radius(IndexType n) const;
    float varianceR(IndexType n) const;
    float varianceZ(IndexType n) const;    

    // The get method in trackcontainer creates every time a new 
    // track state proxy giving this container as input
    // would that be ok for space points? 

  private:
    holder_t<container_t> m_container;
    std::vector<SpacePointProxyType> m_proxies {}; // this will go away ?
  };

  // Deduction rules
  template<typename container_t>
  SpacePointContainer(container_t& container)
    -> SpacePointContainer<container_t, Acts::detail_tc::RefHolder>;

  template<typename container_t>
  SpacePointContainer(container_t&& container)
    -> SpacePointContainer<container_t, Acts::detail_tc::ValueHolder>;
  
  // Implementations
  template<typename container_t, 
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::IndexType
    SpacePointContainer<container_t, holder_t>::size() const
    { return container().size_impl(); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::iterator 
    SpacePointContainer<container_t, holder_t>::begin() 
    { return {*this, 0}; }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::iterator 
    SpacePointContainer<container_t, holder_t>::end() 
    { return {*this, size()}; }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::const_iterator 
    SpacePointContainer<container_t, holder_t>::begin() const 
    { return {*this, 0}; }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::const_iterator 
    SpacePointContainer<container_t, holder_t>::end() const 
    { return {*this, size()}; }
    
  template<typename container_t,
    template <typename> class holder_t>
    inline const container_t& 
    SpacePointContainer<container_t, holder_t>::container() const
    { return *m_container; }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float 
    SpacePointContainer<container_t, holder_t>::x(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const 
    { return container().x_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float SpacePointContainer<container_t, holder_t>::y(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
    { return container().y_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float SpacePointContainer<container_t, holder_t>::z(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
    { return container().z_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float SpacePointContainer<container_t, holder_t>::radius(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
    { return container().radius_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float SpacePointContainer<container_t, holder_t>::varianceR(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
    { return container().varianceR_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline float SpacePointContainer<container_t, holder_t>::varianceZ(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
    { return container().varianceZ_impl(n); }
  
  template<typename container_t,
    template <typename> class holder_t>
    inline typename SpacePointContainer<container_t, holder_t>::SpacePointProxyType&
    SpacePointContainer<container_t, holder_t>::get(
      typename SpacePointContainer<container_t, holder_t>::IndexType n)
  { return m_proxies[n]; }

  template<typename container_t,
    template <typename> class holder_t>
    inline const typename SpacePointContainer<container_t, holder_t>::SpacePointProxyType&
    SpacePointContainer<container_t, holder_t>::get(
      typename SpacePointContainer<container_t, holder_t>::IndexType n) const
  { return m_proxies.at(n); }

}
