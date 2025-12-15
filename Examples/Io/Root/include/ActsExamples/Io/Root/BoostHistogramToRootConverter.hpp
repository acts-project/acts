// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Validation/BoostHistogramWrappers.hpp"

#include <TH1F.h>
#include <TH2F.h>

namespace ActsExamples::BoostHistogramToRoot {

/// Convert BoostHistogram1D to ROOT TH1F
///
/// Creates a new ROOT TH1F histogram with the same binning, content, and
/// errors as the input boost histogram. The histogram uses Sumw2() for proper
/// error tracking with weighted fills.
///
/// @param boostHist The boost histogram to convert
/// @return Raw pointer to new TH1F (caller owns and must delete)
TH1F* toTH1F(const BoostHistogram1D& boostHist);

/// Convert BoostHistogram2D to ROOT TH2F
///
/// Creates a new ROOT TH2F histogram with the same binning, content, and
/// errors as the input boost histogram. The histogram uses Sumw2() for proper
/// error tracking with weighted fills.
///
/// @param boostHist The boost histogram to convert
/// @return Raw pointer to new TH2F (caller owns and must delete)
TH2F* toTH2F(const BoostHistogram2D& boostHist);

}  // namespace ActsExamples::BoostHistogramToRoot
