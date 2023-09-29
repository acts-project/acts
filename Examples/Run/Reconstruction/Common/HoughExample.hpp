// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <memory>

namespace ActsExamples {
class IBaseDetector;
}

/// The Propagation example
///
///
/// @param argc the number of arguments of the call
/// @param argv the argument list
/// @param detector The detector descriptor instance
int runHoughExample(
    int argc, char* argv[],
    const std::shared_ptr<ActsExamples::IBaseDetector>& detector);
