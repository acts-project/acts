// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Result.hpp"
#include "ActsAlignment/Geometry/AlignableStructure.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"

#include <memory>
#include <span>
#include <string>

namespace ActsAlignment::detail {

/// Preparation failure with diagnostic context, without logging an input error.
struct AlignmentTargetError {
  AlignmentError code;
  /// Target/phase and original provider message, when available.
  std::string message;
};

/// Validate the graph, registered placement uniqueness, and explicit selection.
/// Safe before constructing AlignmentHierarchy, including for cyclic input.
/// None selections require no numerical frame and may have
/// children/constraints.
/// @param structures Registered roots
/// @param activeStructures Explicit selections, which must be registered
/// @return Whether any selected structure has active DoFs, or a diagnostic
Acts::Result<bool, AlignmentTargetError> validateAlignmentSelection(
    std::span<const std::shared_ptr<AlignableStructure>> structures,
    std::span<const std::shared_ptr<AlignableStructure>> activeStructures);

}  // namespace ActsAlignment::detail
