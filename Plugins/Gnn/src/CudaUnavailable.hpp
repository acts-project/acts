// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string_view>

namespace ActsPlugins::detail {

/// Report that a CUDA-only code path was reached in a build configured without
/// CUDA. Keeps each `#else` arm of `#ifdef ACTS_GNN_WITH_CUDA` to a single call
/// against one always-compiled definition.
///
/// @param operation what the caller was asked to do, phrased to follow
///        "Cannot ", e.g. "clone CUDA tensor"
[[noreturn]] void throwNoCudaSupport(std::string_view operation);

}  // namespace ActsPlugins::detail
