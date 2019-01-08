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
  struct path_length_sorter
  {
  public:
    template <typename identifier_t, typename parameters_t, typename jacobian_t>
    bool
    operator()(const TrackState<identifier_t, parameters_t, jacobian_t>& lhs,
               const TrackState<identifier_t, parameters_t, jacobian_t>& rhs)
    {
      return lhs.parametric.pathLength < rhs.parametric.pathLength;
    }
  };
}
}
