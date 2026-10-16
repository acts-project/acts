// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Checks that the standard library's own precondition checks are enabled.
// Only libstdc++ is covered: _LIBCPP_DEBUG is a no-op in current libc++.

#include <cstddef>
#include <cstdio>
#include <vector>

int main() {
#if defined(__GLIBCXX__)
  std::vector<int> values(1);
  volatile std::size_t index = 1;
  int value = values[index];
  std::printf("ACTS_ASSERTION_NOT_TRIGGERED %d\n", value);
  return 0;
#else
  // Reported as skipped by RunAssertProbe.cmake
  return 77;
#endif
}
