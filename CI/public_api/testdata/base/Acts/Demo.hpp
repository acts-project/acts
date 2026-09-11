// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file fixture header for the API-surface self-test (base revision)

namespace Acts {

/// A demo type.
struct Demo {
  /// A tolerance.
  double tolerance = 0.0;
  /// A label.
  int label = 0;

  /// Do a thing.
  /// @param x the x
  void doThing(int x);

  /// A function that will be removed.
  /// @param a the a
  void oldFn(double a);
};

/// A free function.
/// @param n a number
int compute(int n);

}  // namespace Acts
