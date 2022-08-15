// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Set the Calibration Context PLUGIN
#ifdef ACTS_CORE_CALIBRATIONCONTEXT_PLUGIN
#include ACTS_CORE_CALIBRATIONCONTEXT_PLUGIN
#else

#include <any>

namespace Acts {

/// @brief This is the central definition of the Acts
/// payload object regarding detector calbiration
///
/// It is propagated through the code to allow for event/thread
/// dependent calibration

using CalibrationContext = std::any;

}  // namespace Acts

#endif
