// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Histogram.hpp"

#include <vector>

class TEfficiency;
class TH1F;
class TH2F;
class TProfile;

namespace ActsPlugins {

/// Convert Histogram1 to ROOT TH1F
///
/// Creates a new ROOT TH1F histogram with the same binning, content, and
/// errors as the input boost histogram.
///
/// @param boostHist The boost 1D histogram to convert
/// @return Raw pointer to new TH1F (caller owns and must delete)
TH1F* toRoot(const Acts::Experimental::Histogram1& boostHist);

/// Convert Histogram2 to ROOT TH2F
///
/// Creates a new ROOT TH2F histogram with the same binning, content, and
/// errors as the input boost histogram.
///
/// @param boostHist The boost 2D histogram to convert
/// @return Raw pointer to new TH2F (caller owns and must delete)
TH2F* toRoot(const Acts::Experimental::Histogram2& boostHist);

/// Convert ProfileHistogram1 to ROOT TProfile
///
/// Creates a new ROOT TProfile histogram with the same binning and mean values
/// as the input boost profile histogram.
///
/// @param boostProfile The boost profile histogram to convert
/// @return Raw pointer to new TProfile (caller owns and must delete)
TProfile* toRoot(const Acts::Experimental::ProfileHistogram1& boostProfile);

/// Convert Efficiency1 to ROOT TEfficiency
///
/// Creates a new ROOT TEfficiency object with the same binning, passed counts,
/// and total counts as the input boost efficiency histogram.
///
/// @param boostEff The boost 1D efficiency histogram to convert
/// @return Raw pointer to new TEfficiency (caller owns and must delete)
TEfficiency* toRoot(const Acts::Experimental::Efficiency1& boostEff);

/// Convert Efficiency2 to ROOT TEfficiency
///
/// Creates a new ROOT TEfficiency object (2D) with the same binning, passed
/// counts, and total counts as the input boost efficiency histogram.
///
/// @param boostEff The boost 2D efficiency histogram to convert
/// @return Raw pointer to new TEfficiency (caller owns and must delete)
TEfficiency* toRoot(const Acts::Experimental::Efficiency2& boostEff);

}  // namespace ActsPlugins
