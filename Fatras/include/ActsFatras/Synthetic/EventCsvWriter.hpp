// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <string>

namespace ActsFatras::Synthetic {

/// Write an event to two CSV files, `<prefix>_spacepoints.csv` and
/// `<prefix>_particles.csv`, for offline inspection of the distributions.
///
/// @param event the generated event
/// @param layout the detector the event was generated on
/// @param prefix the path prefix of the two files
void writeEventCsv(const Event& event, const DetectorLayout& layout,
                   const std::string& prefix);

}  // namespace ActsFatras::Synthetic
