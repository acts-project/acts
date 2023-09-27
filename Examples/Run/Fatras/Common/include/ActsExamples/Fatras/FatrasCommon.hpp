// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

namespace ActsExamples {
class IBaseDetector;
}

/// Main function for running Fatras with a specific detector.
///
/// @param argc number of command line arguments
/// @param argv command line arguments
/// @param detector is the detector to be used
int runFatras(int argc, char* argv[],
              const std::shared_ptr<ActsExamples::IBaseDetector>& detector);
