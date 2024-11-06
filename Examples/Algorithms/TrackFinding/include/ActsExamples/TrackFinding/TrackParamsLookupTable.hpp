// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

namespace ActsExamples {

/// @brief Track parameters lookup table axis used
/// in the track estimation algorithm
using LookupAxis =
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;

/// @brief Track parameters lookup table axis generator
/// used in the track estimation algorithm
using LookupAxisGen = Acts::GridAxisGenerators::EqOpenEqOpen;

/// @brief Grid used to accumulate IP track parameters and
/// reference layer track parameters for a given position
/// in the track estimation algorithm
using LookupAccumGrid =
    Acts::Grid<std::vector<Acts::CurvilinearTrackParameters>, LookupAxis,
               LookupAxis>;

/// @brief Container for the lookup accumulation grids to
/// handle multiple reference layers in the track estimation
using LookupAccumGridContainer =
    std::unordered_map<Acts::GeometryIdentifier, LookupAccumGrid>;

/// @brief IP-reference layer track parameters lookup pair
using LookupPair = std::pair<std::shared_ptr<Acts::CurvilinearTrackParameters>,
                             std::shared_ptr<Acts::CurvilinearTrackParameters>>;

/// @brief Lookup grid for track parameters estimation
/// in a given layer
using LookupGrid = Acts::Grid<LookupPair, LookupAxis, LookupAxis>;

/// @brief Lookup table for track parameters estimation
/// in the track estimation algorithm
using Lookup = std::unordered_map<Acts::GeometryIdentifier, LookupGrid>;

}  // namespace ActsExamples
