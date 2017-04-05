// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <tuple>
#include "ACTS/Utilities/detail/get_nbins_helper.hpp"
#include "ACTS/Utilities/detail/global_bin_helper.hpp"

namespace Acts {

namespace Test {
  template <typename T, class... Axes>
  struct GridTester;
}

namespace detail {

  template <typename T, class... Axes>
  class Grid
  {
    /// unit test helper class
    friend struct Test::GridTester<T, Axes...>;

    static constexpr size_t dimension = sizeof...(Axes);

  public:
    typedef T                 value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;

    Grid(std::tuple<Axes...> axes) : m_axes(std::move(axes))
    {
      m_values.resize(size());
    }

    template <class Point>
    reference
    at(const Point& point)
    {
      return m_values.at(getGlobalBinIndex(point));
    }

    template <class Point>
    const_reference
    at(const Point& point) const
    {
      return m_values.at(getGlobalBinIndex(point));
    }

    size_t
    size() const
    {
      return get_nbins_helper::getNBins(m_axes);
    }

  private:
    template <class Point>
    size_t
    getGlobalBinIndex(const Point& point) const
    {
      return global_bin_helper::getGlobalBin(point, m_axes);
    }

    std::array<size_t, dimension>
    getLocalBinIndices(size_t bin) const
    {
      return global_bin_helper::getLocalBinIndices(bin, m_axes);
    }

    std::tuple<Axes...> m_axes;
    std::vector<T>      m_values;
  };
}  // namespace detail

}  // namespace Acts
