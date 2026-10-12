// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "hip_error_handling.hpp"

// System include(s).
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace traccc::hip::details {

void throw_error(hipError_t errorCode, const char* expression, const char* file,
                 int line) {
  // Create a nice error message.
  std::ostringstream errorMsg;
  errorMsg << file << ":" << line << " Failed to execute: " << expression
           << " (" << hipGetErrorString(errorCode) << ")";

  // Now throw a runtime error with this message.
  throw std::runtime_error(errorMsg.str());
}

}  // namespace traccc::hip::details
