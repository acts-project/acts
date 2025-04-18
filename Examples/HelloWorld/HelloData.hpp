// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <ostream>
#include <vector>

namespace ActsExamples {

/// An example data object to be shared via the event store.
struct HelloData {
  double x;
  double a;
  double b;
  double t;
};

inline std::ostream& operator<<(std::ostream& os, const HelloData& data) {
  os << data.x << ", " << data.a << ", " << data.b << ", " << data.t;
  return os;
}

using HelloDataCollection = std::vector<HelloData>;

}  // namespace ActsExamples
