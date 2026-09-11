// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// @file fixture header for the API-surface self-test (head revision)

namespace Acts {

/// A demo type.
struct Demo {
  /// A tolerance (now single precision -- BREAKING retype).
  float tolerance = 0.0;
  /// A label.
  int label = 0;
  /// A newly added member (addition, not breaking).
  bool flag = false;

  /// Do a thing, now with an extra required argument (BREAKING).
  /// @param x the x
  /// @param y the y
  void doThing(int x, int y);

  // oldFn(double) removed -- BREAKING.
};

/// A free function, now with a *defaulted* argument (addition, not breaking).
/// @param n a number
/// @param scale a scale factor
int compute(int n, double scale = 1.0);

/// A brand new type (addition).
struct NewThing {};

}  // namespace Acts
