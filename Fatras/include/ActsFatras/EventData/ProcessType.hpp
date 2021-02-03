// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>
#include <iosfwd>

namespace ActsFatras {

/// Process type identifier.
///
/// Encodes the type of process that generated a particle.
enum class ProcessType : uint32_t {
  eUndefined = 0,
};

std::ostream &operator<<(std::ostream &os, ProcessType processType);

}  // namespace ActsFatras
