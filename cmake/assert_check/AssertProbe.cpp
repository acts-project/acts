// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Built with NDEBUG defined: the assert below only fires if
// ACTS_FORCE_ASSERTIONS re-enabled it.

#include <cassert>
#include <cstdio>

int main() {
  assert(false && "forced assertion");
  std::puts("ACTS_ASSERTION_NOT_TRIGGERED");
  return 0;
}
