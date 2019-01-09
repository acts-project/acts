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
/// @brief  Sorter for boost_variant
struct TrackStatePathLengthSorter
{
public:
  template <typename identifier_t, typename parameters_t>
  bool
  operator()(const TrackState<identifier_t, parameters_t>& lhs,
             const TrackState<identifier_t, parameters_t>& rhs)
  {
    return lhs.parameter.pathLength < rhs.parameter.pathLength;
  }
};
}
