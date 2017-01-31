// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <tuple>
#include <utility>

namespace Acts {

namespace detail {

  template <size_t N>
  struct get_nbins_helper_impl;

  template <size_t N>
  struct get_nbins_helper_impl
  {
    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      size_t thisAxisNBins = std::get<N>(axes).getNBins();
      return thisAxisNBins * get_nbins_helper_impl<N - 1>::getNBins(axes);
    }
  };

  template <>
  struct get_nbins_helper_impl<0u>
  {
    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      size_t thisAxisNBins = std::get<0u>(axes).getNBins();
      return thisAxisNBins;
    }
  };

  struct get_nbins_helper
  {
    template <class... Axes>
    static size_t
    getNBins(const std::tuple<Axes...>& axes)
    {
      constexpr size_t MAX = sizeof...(Axes)-1;
      return get_nbins_helper_impl<MAX>::getNBins(axes);
    }
  };

}  // namespace detail

}  // namespace Acts
