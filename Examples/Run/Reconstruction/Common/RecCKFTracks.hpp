// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>

namespace ActsExamples {
class IBaseDetector;
}

/// The options for running CKF reco
///
/// @param desc The options description to add options to
void addRecCKFOptions(ActsExamples::Options::Description& desc);

/// Main function for running CKF reco with a specific detector.
///
/// @param argc number of command line arguments
/// @param argv command line arguments
/// @param detector is the detector to be used
int runRecCKFTracks(
    int argc, char* argv[],
    const std::shared_ptr<ActsExamples::IBaseDetector>& detector);
