// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// boost include(s)
#include <boost/variant.hpp>

namespace Acts {
namespace detail {
  /// @brief  Sorter for boost_variant
  struct variant_path_length_sorter : public boost::static_visitor<bool>
  {
  public:
    /// @brief compares different types
    template <typename type_a_t, typename type_b_t>
    bool
    operator()(const type_a_t& lhs, const type_b_t& rhs) const
    {
      return (lhs.parametric.pathLength < rhs.parametric.pathLength);
    }
  };

  /// @brief vistor pattern path length sorter
  struct path_length_sorter
  {
  public:
    template <typename variant_type_t>
    bool
    operator()(const variant_type_t& lhs, const variant_type_t& rhs)
    {
      return boost::apply_visitor(variant_path_length_sorter(), lhs, rhs);
    }
  };
}
}