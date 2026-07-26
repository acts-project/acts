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
