// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "CudaUnavailable.hpp"

#include <stdexcept>
#include <string>

namespace ActsPlugins::detail {

void throwNoCudaSupport(std::string_view operation) {
  throw std::runtime_error("Cannot " + std::string(operation) +
                           ": ACTS was not compiled with CUDA support");
}

}  // namespace ActsPlugins::detail
