// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts::detail {
/// Enum which determines how the axis handle its outer boundaries
/// possible values values
/// - Open is the default behaviour: out of bounds
/// positions are filled into the over or underflow bins
/// - Bound: out-of-bounds positions resolve to first/last bin
/// respectively
/// - Closed: out-of-bounds positions resolve to the outermost
/// bin on the opposite side
enum class AxisBoundaryType { Open, Bound, Closed };

/// Enum which determines the binning type of the axis
enum class AxisType { Equidistant, Variable };

/// @brief calculate bin indices from a given binning structure
///
/// This class provides some basic functionality for calculating bin indices
/// for a given binning configuration. Both equidistant as well as variable
/// binning structures are supported.
///
/// Bin intervals are defined such that the lower bound is closed and the
/// upper bound is open.
///
/// @tparam equidistant flag whether binning is equidistant (@c true)
///                     or not (@c false)
template <AxisType type, AxisBoundaryType bdt = AxisBoundaryType::Open>
class Axis;

using EquidistantAxis = Axis<AxisType::Equidistant>;
using VariableAxis = Axis<AxisType::Variable>;

}  // namespace Acts::detail
