// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"

#include <nlohmann/json.hpp>

namespace ActsPlugins {

/// Convert Histogram1 to JSON
///
/// Schema: {"name", "title", "type":"histogram", "axes":[...], "values":[...]}
/// Axes carry {"edges":[...], "label":"..."}.
/// Values are inner bins only (no flow), in row-major order.
///
/// @param boostHist The 1D histogram to convert
/// @return JSON representation
nlohmann::json toJson(const Acts::Experimental::Histogram1& boostHist);

/// Convert Histogram2 to JSON
///
/// @param boostHist The 2D histogram to convert
/// @return JSON representation
nlohmann::json toJson(const Acts::Experimental::Histogram2& boostHist);

/// Convert Histogram3 to JSON
///
/// @param boostHist The 3D histogram to convert
/// @return JSON representation
nlohmann::json toJson(const Acts::Experimental::Histogram3& boostHist);

/// Convert ProfileHistogram1 to JSON
///
/// Schema: {"name", "title", "type":"profile", "axes":[...],
///          "counts":[...], "means":[...], "sum_of_deltas_squared":[...]}
/// sampleAxisTitle is stored under "sampleAxisTitle".
///
/// @param boostProfile The 1D profile histogram to convert
/// @return JSON representation
nlohmann::json toJson(
    const Acts::Experimental::ProfileHistogram1& boostProfile);

/// Convert Efficiency1 to JSON
///
/// Schema: {"name", "title", "type":"efficiency", "axes":[...],
///          "accepted":[...], "total":[...]}
///
/// @param boostEff The 1D efficiency histogram to convert
/// @return JSON representation
nlohmann::json toJson(const Acts::Experimental::Efficiency1& boostEff);

/// Convert Efficiency2 to JSON
///
/// @param boostEff The 2D efficiency histogram to convert
/// @return JSON representation
nlohmann::json toJson(const Acts::Experimental::Efficiency2& boostEff);

}  // namespace ActsPlugins
