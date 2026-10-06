// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace traccc {
template <typename T1, typename T2>
struct pair {
 public:
  using first_type = T1;
  using second_type = T2;

  T1 first;
  T2 second;
};
}  // namespace traccc
