// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/Holders.hpp"

#include "ActsExamples/EventData/SimSpacePoint.hpp"

// memory backend is a vector of pointer of SimSpacePoint

namespace ActsExamples {
  
  template<typename collection_t>
  class SpacePointContainer {
  public:
    using CollectionType = collection_t;
    using ValueType = typename CollectionType::value_type;

    friend Acts::SpacePointContainer<ActsExamples::SpacePointContainer<collection_t>, Acts::detail::RefHolder>;
    
    // default constructor is of no use. It cannot be used, so why bother?
    SpacePointContainer() = delete;
    // we never get the ownership. In both read-only and read-and-write mode
    // the memory backend is independetly handled. This is only interfacing it to ACTS
    SpacePointContainer(CollectionType&& container) = delete;
    SpacePointContainer(CollectionType& container)
      : m_storage(container)
    {}
    SpacePointContainer(CollectionType* container)
      : m_storage(container)
    {}
    
    // No copy contructor or copy operation allowed
    SpacePointContainer(const SpacePointContainer<collection_t>&) = delete;
    SpacePointContainer<collection_t>& operator=(const SpacePointContainer<collection_t>&) = delete;

    // only move operation allowed
    SpacePointContainer(SpacePointContainer<collection_t>&& other) noexcept
      : m_storage( std::exchange(other.m_storage.ptr, nullptr) )
    {}
    SpacePointContainer<collection_t>& operator=(SpacePointContainer<collection_t>&& other) noexcept
    {
      m_storage = std::exchange(other.m_storage.ptr, nullptr);
      return *this;
    }
    
    ~SpacePointContainer() = default;

  private:
    using IndexType = std::size_t;

    IndexType size_impl() const;
    float x_impl(IndexType idx) const;
    float y_impl(IndexType idx) const;
    float z_impl(IndexType idx) const;
    float radius_impl(IndexType idx) const;
    float varianceR_impl(IndexType idx) const;
    float varianceZ_impl(IndexType idx) const;

    // template<typename T>
    // const T& component(HashedString key, std::size_t n) const
    // {
    //   using namespace Acts::HashedStringLiteral;
    //   switch (key) {
    //   case ""_hash:
    // 	return *std::any_cast<T*>(container().component(HashedString(key), n));
    //   default:
    // 	throw std::runtime_error("no such component " + std::to_string(key));
    //   }
    // }
    
  private:
    const CollectionType& storage() const;
    
  private:
    Acts::detail::RefHolder<CollectionType> m_storage;
  };

  template<typename collection_t>  
  inline typename SpacePointContainer<collection_t>::IndexType
  SpacePointContainer<collection_t>::size_impl() const
  { return storage().size(); }

  // TO-DO
  // Be smart here... collection_t can container values or pointers ...
  
  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::x_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  { return storage()[idx]->x(); }

  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::y_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  { return storage()[idx]->y(); }

  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::z_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  { return storage()[idx]->z(); }

  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::radius_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  {
    const float x = x_impl(idx);
    const float y = y_impl(idx);
    return std::sqrt(x*x + y*y);
  }

  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::varianceR_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  { return storage()[idx]->varianceR(); }

  template<typename collection_t>
  inline float
  SpacePointContainer<collection_t>::varianceZ_impl(typename SpacePointContainer<collection_t>::IndexType idx) const
  { return storage()[idx]->varianceZ(); }

  template<typename collection_t>
  const typename SpacePointContainer<collection_t>::CollectionType&
  SpacePointContainer<collection_t>::storage() const
  { return *m_storage; }

}
