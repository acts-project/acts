// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace traccc {

// Simple struct for combining multiple memory resources
struct memory_resource {
  // device or shared memory resource
  vecmem::memory_resource& main;

  // optional host accessible memory resource
  vecmem::memory_resource* host = nullptr;
};

}  // namespace traccc
