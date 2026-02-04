// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonSegment.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <string>
#include <vector>

namespace ActsExamples {

/// @brief Visualizes a Hough transform accumulator plane and found maxima
///
/// Creates a PDF showing:
/// - The Hough accumulator plane as a 2D histogram
/// - Found maxima with uncertainty boxes (blue)
/// - Truth segment parameters (red crosses)
///
/// @param outputPath Full path for the output PDF file
/// @param bucketId Identifier of the station bucket being visualized
/// @param maxima Found Hough maxima from the peak finder
/// @param plane Filled Hough accumulator plane
/// @param axis Axis ranges of the Hough plane
/// @param truthSegments Container of truth segments for comparison
/// @param logger Logger for diagnostic output
void visualizeMuonHoughMaxima(
    const std::string& outputPath, const MuonSpacePoint::MuonId& bucketId,
    const std::vector<Acts::HoughTransformUtils::PeakFinders::
                          IslandsAroundMax<const MuonSpacePoint*>::Maximum>&
        maxima,
    const Acts::HoughTransformUtils::HoughPlane<const MuonSpacePoint*>& plane,
    const Acts::HoughTransformUtils::HoughAxisRanges& axis,
    const MuonSegmentContainer& truthSegments, const Acts::Logger& logger);

}  // namespace ActsExamples
