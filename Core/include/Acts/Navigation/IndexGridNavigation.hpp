// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"

#include <algorithm>
#include <array>

namespace Acts {

/// @brief  This is an index grid based navigation state updator, it uses
/// an extractor type and a filler type to handle the navigation state
///
/// @note a transform is applied `p3l = transform * p3` in order to allow
/// shifted, transformed grids
///
/// It can be used for volumes, surfaces at convenience
///
/// @tparam navigation_type distinguishes between internal and external navigation
/// @tparam grid_t is the type of the grid
/// @tparam extractor_type is the helper to extract the object
/// @tparam filler_type is the helper to fill the object into the nState
template <typename grid_t>
class IndexGridNavigation {
 public:
  /// Broadcast the grid type
  using grid_type = grid_t;

  /// The grid where the indices are stored
  grid_type grid;

  /// These are the cast parameters - copied from constructor
  std::array<AxisDirection, grid_type::DIM> casts{};

  /// A transform to be applied to the position
  Transform3 transform = Transform3::Identity();

  /// @brief  Constructor for a grid based surface attacher
  /// @param igrid the grid that is moved into this attacher
  /// @param icasts is the cast values array
  /// @param itr a transform applied to the global position
  IndexGridNavigation(grid_type&& igrid,
                      const std::array<AxisDirection, grid_type::DIM>& icasts,
                      const Transform3& itr = Transform3::Identity())
      : grid(std::move(igrid)), casts(icasts), transform(itr) {}

  IndexGridNavigation() = delete;
};

}  // namespace Acts
