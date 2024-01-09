// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsExamples {
class IBaseDetector;
}

/// The Propagation example
///
///
/// @param argc the number of arguments of the call
/// @param argv the argument list
/// @param detector The detector descriptor instance
int propagationExample(int argc, char* argv[],
                       ActsExamples::IBaseDetector& detector);
