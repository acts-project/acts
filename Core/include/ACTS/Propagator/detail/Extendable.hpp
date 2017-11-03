// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EXTENDABLE_HPP
#define ACTS_EXTENDABLE_HPP 1

#include <tuple>
#include <type_traits>
#include "ACTS/Utilities/detail/MPL/all_of.hpp"

namespace Acts {

namespace detail {

  template <typename... Extensions>
  struct Extendable : private std::tuple<Extensions...>
  {
    // clang-format off
    static_assert(detail::all_of_v<std::is_default_constructible<Extensions>::value...>,
                  "all extensions must be default constructible");
    static_assert(detail::all_of_v<std::is_copy_constructible<Extensions>::value...>,
                  "all extensions must be copy constructible");
    static_assert(detail::all_of_v<std::is_move_constructible<Extensions>::value...>,
                  "all extensions must be move constructible");
    // clang-format on

    // Explicit constructor
    explicit Extendable(const Extensions&... extensions)
      : std::tuple<Extensions...>(extensions...)
    {
    }

    template <typename R>
    const R&
    get() const
    {
      return std::get<R>(*this);
    }

    template <typename R>
    R&
    get()
    {
      return std::get<R>(*this);
    }

    const std::tuple<Extensions...>&
    tuple() const
    {
      return (*this);
    }

    std::tuple<Extensions...>&
    tuple()
    {
      return (*this);
    }
  };

}  // namespace detail

}  // namespace Acts
#endif  // ACTS_EXTENDABLE_HPP
