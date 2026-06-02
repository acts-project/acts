// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

namespace ActsPython {

std::unordered_map<PyObject*, WhiteBoardRegistry::RegistryEntry>&
WhiteBoardRegistry::instance() {
  static std::unordered_map<PyObject*, RegistryEntry> map;
  return map;
}

}  // namespace ActsPython
