// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <memory>
#include <unordered_map>

namespace ActsExamples {

using TrackParamsLookupPair =
    std::pair<std::shared_ptr<Acts::CurvilinearTrackParameters>,
              std::shared_ptr<Acts::CurvilinearTrackParameters>>;

/// @brief Track parameters lookup table axis used
/// in the track estimation algorithm
using TrackParamsLookupAxis =
    Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;

/// @brief Track parameters lookup table axis generator
/// used in the track estimation algorithm
using TrackParamsLookupAxisGen = Acts::GridAxisGenerators::EqOpenEqOpen;

/// @brief Lookup grid for track parameters estimation
/// in a given layer
using TrackParamsLookupGrid =
    Acts::Grid<TrackParamsLookupPair, TrackParamsLookupAxis,
               TrackParamsLookupAxis>;

/// @brief Lookup table for track parameters estimation
/// in the track estimation algorithm
using TrackParamsLookup =
    std::unordered_map<Acts::GeometryIdentifier, TrackParamsLookupGrid>;

}  // namespace ActsExamples
