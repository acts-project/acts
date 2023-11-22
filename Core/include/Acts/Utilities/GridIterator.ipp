// -*- C++ -*-
// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>

namespace Acts {
  // Global Iterator
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>::GridGlobalIterator(const Acts::Grid<T, Axes ...>& grid,
						      std::size_t idx)
    :  m_grid( &grid ),
       m_idx( idx )
  {}
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>::GridGlobalIterator(GridGlobalIterator<T, Axes ...>&& other) noexcept
    : m_grid( std::exchange(other.m_grid.ptr, nullptr) ),
      m_idx( other.m_idx )
  {}
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>&
  GridGlobalIterator<T, Axes ...>::operator=(GridGlobalIterator<T, Axes ...>&& other) noexcept
  {
    m_grid.ptr = std::exchange(other.m_grid.ptr, nullptr);
    m_idx = other.m_idx;
    return *this;
  }

  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator==(const GridGlobalIterator<T, Axes ...>& other) const
  {
    return (m_grid.ptr == other.m_grid.ptr) && m_idx == other.m_idx;
  }
  
  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator!=(const GridGlobalIterator<T, Axes ...>& other) const
  {
    return !(*this == other);
  }
  
  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator<(const GridGlobalIterator<T, Axes ...>& other) const
  { return m_idx < other.m_idx; }
  
  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator>(const GridGlobalIterator<T, Axes ...>& other) const
  { return m_idx > other.m_idx;  }
  
  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator<=(const GridGlobalIterator<T, Axes ...>& other) const
  { return !(*this > other); }
  
  template <typename T, class ... Axes>
  bool GridGlobalIterator<T, Axes ...>::operator>=(const GridGlobalIterator<T, Axes ...>& other) const
  { return !(*this < other); }
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>&
  GridGlobalIterator<T, Axes ...>::operator+=(const std::size_t offset)
  {
    m_idx += offset;
    return *this;
  }
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>&
  GridGlobalIterator<T, Axes ...>::operator-=(const std::size_t offset)
  {
    m_idx -= offset;
    return *this;
  }
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>
  GridGlobalIterator<T, Axes ...>::operator+(const std::size_t offset) const
  {
    return {m_grid.ptr, m_idx + offset};
  }
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>
  GridGlobalIterator<T, Axes ...>::operator-(const std::size_t offset) const
  {
    return {m_grid.ptr, m_idx - offset};
  }
  
  template <typename T, class ... Axes>
  typename GridGlobalIterator<T, Axes ...>::difference_type
  GridGlobalIterator<T, Axes ...>::operator-(const GridGlobalIterator<T, Axes ...>& other) const {
    assert(other > *this);
    return other.m_idx - m_idx;
  }
  
  template <typename T, class ... Axes>
  const typename GridGlobalIterator<T, Axes ...>::value_type&
  GridGlobalIterator<T, Axes ...>::operator*() const
  {
    return m_grid->at(m_idx);
  }
  
  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>&
  GridGlobalIterator<T, Axes ...>::operator++() {
    ++m_idx;
    return *this;
  }

  template <typename T, class ... Axes>
  GridGlobalIterator<T, Axes ...>
  GridGlobalIterator<T, Axes ...>::operator++(int) {
    GridGlobalIterator<T, Axes ...> output(m_grid, m_idx++);    
    return output;
  }
  

  

     

  // Local Iterator
  template <typename T, class... Axes>
  Acts::GridLocalIterator<T, Axes ...>::GridLocalIterator(const Acts::Grid<T, Axes ...>& grid,
							  std::array<std::size_t, DIM> indexes)
    : m_grid( &grid ),
      m_numLocalBins( std::move(grid.numLocalBins()) ),
      m_currentIndex( std::move(indexes) )
  {
    for (std::size_t i(0); i<DIM; ++i) {
      m_navigationIndex[i].resize(m_numLocalBins[i]);
      std::iota(m_navigationIndex[i].begin(), m_navigationIndex[i].end(), 1ul);
      m_localPosition[i] = 1ul;
    }
  }

  template <typename T, class... Axes>
   Acts::GridLocalIterator<T, Axes ...>::GridLocalIterator(const Acts::Grid<T, Axes ...>& grid,
							   std::array<std::size_t, DIM> indexes,
							   std::array<std::vector<std::size_t>, DIM> navigation)
    : m_grid( &grid ),
      m_numLocalBins( std::move(grid.numLocalBins()) ),
      m_currentIndex( std::move(indexes) ),
      m_navigationIndex( std::move(navigation) )
  {
    // check navigation consistency
    for (std::size_t i(0); i<DIM; ++i) {
      if (m_navigationIndex[i].size() != m_numLocalBins[i]) {
	throw std::invalid_argument("Invalid navigation sequence in local grid iterator.");
      }
      m_localPosition[i] = m_navigationIndex[i][m_currentIndex[i]];
    }
  }
  
  template <typename T, class... Axes>
  Acts::GridLocalIterator<T, Axes ...>::GridLocalIterator(Acts::GridLocalIterator<T, Axes ...>&& other) noexcept
  : m_grid( std::exchange(other.m_grid.ptr, nullptr) ),
    m_numLocalBins( std::move(other.m_numLocalBins) ),
    m_currentIndex( std::move(other.m_currentIndex) ),
    m_navigationIndex( std::move(other.m_navigationIndex) ),
    m_localPosition( std::move(other.m_localPosition) )
  {}
  
  template <typename T, class... Axes>
  Acts::GridLocalIterator<T, Axes ...>&
  Acts::GridLocalIterator<T, Axes ...>::operator=(Acts::GridLocalIterator<T, Axes ...>&& other) noexcept
  {
    m_grid.ptr = std::exchange(other.m_grid.ptr, nullptr);
    m_numLocalBins = std::move(other.m_numLocalBins);
    m_currentIndex = std::move(other.m_currentIndex);
    m_navigationIndex = std::move(other.m_navigationIndex);
    m_localPosition = std::move(other.m_localPosition);
    return *this;
  }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator==(const Acts::GridLocalIterator<T, Axes ...>& other) const
  {
    if (m_grid.ptr != other.m_grid.ptr) {
      return false;
    }
    
    for (std::size_t i(0); i<DIM; ++i) {
      if (m_currentIndex[i] != other.m_currentIndex[i]) {
	return false;
      }
    }
    
    return true;
  }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator!=(const Acts::GridLocalIterator<T, Axes ...>& other) const
  { return ! (*this == other); }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator<(const GridLocalIterator<T, Axes ...>& other) const
  { return m_grid.globalBinFromLocalBins(m_currentIndex) < other.m_grid.globalBinFromLocalBins(m_currentIndex); }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator>(const GridLocalIterator<T, Axes ...>& other) const
  { return m_grid.globalBinFromLocalBins(m_currentIndex) > other.m_grid.globalBinFromLocalBins(m_currentIndex); }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator<=(const GridLocalIterator<T, Axes ...>& other) const
  { return ! (*this > other); }
  
  template <typename T, class... Axes>
  bool Acts::GridLocalIterator<T, Axes ...>::operator>=(const GridLocalIterator<T, Axes ...>& other) const
  { return ! (*this < other); }
  
  template <typename T, class... Axes>
  typename Acts::GridLocalIterator<T, Axes ...>::difference_type
  Acts::GridLocalIterator<T, Axes ...>::operator-(const Acts::GridLocalIterator<T, Axes ...>& other) const
  { return other.m_grid->globalBinFromLocalBins(m_currentIndex) - m_grid->globalBinFromLocalBins(m_currentIndex); }

  template <typename T, class... Axes>
  const typename Acts::GridLocalIterator<T, Axes ...>::value_type&
  Acts::GridLocalIterator<T, Axes ...>::operator*() const
  {
    std::array<std::size_t, DIM> localPositionBin;
    for (std::size_t i(0); i<DIM; ++i) {
      localPositionBin[i] = m_navigationIndex[i][m_currentIndex[i]];
    }
    return m_grid->atLocalBins(localPositionBin);
  }

  template <typename T, class... Axes>
  GridLocalIterator<T, Axes ...>&
  GridLocalIterator<T, Axes ...>::operator++()
  {
    increment<DIM - 1>();
    return *this;
  }

  template <typename T, class... Axes>
   GridLocalIterator<T, Axes ...>
  GridLocalIterator<T, Axes ...>::operator++(int)
  {
    GridLocalIterator<T, Axes ...> output(m_grid.ptr, m_currentIndex);
    this->operator++();
    return output;
  }

  template <typename T, class... Axes>
  template <std::size_t N>
  void GridLocalIterator<T, Axes ...>::increment() {
    if (++m_currentIndex[N] < m_numLocalBins[N]) return;
    m_currentIndex[N] = 0;
    if constexpr (N != 0) {
      increment<N-1>();
    } else {
      m_currentIndex = m_numLocalBins;
    }
  }

  template <typename T, class... Axes>
  std::array<std::size_t, GridLocalIterator<T, Axes ...>::DIM>
  GridLocalIterator<T, Axes ...>::localPosition() const
  {
    std::array<std::size_t, DIM> output;
    for (std::size_t i(0); i<DIM; ++i) {
      output[i] = m_navigationIndex[i][m_currentIndex[i]];
    }
    return output;
  }
  
} // namespace Acts

