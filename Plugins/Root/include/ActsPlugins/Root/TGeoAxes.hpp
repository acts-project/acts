// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <ostream>
#include <stdexcept>
#include <string_view>

namespace ActsPlugins {

/// @addtogroup root_plugin
/// @{

/// Type-safe representation of TGeo surface axes.
///
/// The axes string encodes which TGeo coordinate axes map to the local X and Y
/// directions of an Acts surface, with case controlling the sign:
///   - Uppercase (X, Y, Z): positive direction
///   - Lowercase (x, y, z): negative direction
///
/// Valid values are all 3-character permutations of {X,x}, {Y,y}, {Z,z} with
/// no axis letter repeated (e.g. "XYZ", "XZY", "xYz").
///
/// String literals are validated at compile time via the @c consteval
/// constructor — invalid characters or repeated axes are a compile error.
/// Runtime construction is available via @c parse().
class TGeoAxes {
 public:
  /// Construct from a 3-character string literal.
  ///
  /// The parameter type @c const char (&)[4] only binds to exactly 3-character
  /// string literals; wrong lengths are a type error. Invalid characters or
  /// repeated axis letters cause a compile error.
  ///
  /// @param s A 3-character axes string literal, e.g. @c "XYZ"
  consteval TGeoAxes(const char (&s)[4])  // NOLINT(google-explicit-constructor)
      : TGeoAxes(makeChecked(s[0], s[1], s[2])) {}

  /// Construct from a runtime string. Throws on invalid input.
  ///
  /// @param s A string of exactly 3 characters from [XxYyZz] with no repeated
  ///          axis letter.
  /// @throws std::invalid_argument if the string is not a valid axes encoding.
  /// @return The parsed axes as a TGeoAxes object.
  static TGeoAxes parse(std::string_view s) {
    if (s.size() != 3 || !isValid(s[0], s[1], s[2])) {
      throw std::invalid_argument(
          "TGeoAxes: expected exactly 3 characters from [XxYyZz] "
          "with each axis (X, Y, Z) appearing exactly once");
    }
    return TGeoAxes{Checked{s[0], s[1], s[2]}};
  }

  /// Return the axes string as a 3-character string view.
  /// @return The axes string as a 3-character string view.
  std::string_view value() const { return {m_axes, 3}; }

  friend std::ostream& operator<<(std::ostream& os, const TGeoAxes& axes) {
    return os << axes.value();
  }

 private:
  struct Checked {
    char a, b, c;
  };

  static constexpr bool isAxisChar(char c) {
    return c == 'X' || c == 'x' || c == 'Y' || c == 'y' || c == 'Z' || c == 'z';
  }

  static constexpr char baseAxis(char c) {
    if (c == 'x') {
      return 'X';
    }
    if (c == 'y') {
      return 'Y';
    }
    if (c == 'z') {
      return 'Z';
    }
    return c;
  }

  // Shared predicate — constexpr so usable at both compile and runtime.
  static constexpr bool isValid(char a, char b, char c) {
    return isAxisChar(a) && isAxisChar(b) && isAxisChar(c) &&
           baseAxis(a) != baseAxis(b) && baseAxis(b) != baseAxis(c) &&
           baseAxis(a) != baseAxis(c);
  }

  // consteval: a throw reached during evaluation makes the call ill-formed,
  // which manifests as a compile error at the call site.
  static consteval Checked makeChecked(char a, char b, char c) {
    if (!isAxisChar(a) || !isAxisChar(b) || !isAxisChar(c)) {
      throw "TGeoAxes: each character must be one of X x Y y Z z";
    }
    if (baseAxis(a) == baseAxis(b) || baseAxis(b) == baseAxis(c) ||
        baseAxis(a) == baseAxis(c)) {
      throw "TGeoAxes: each axis (X, Y, Z) must appear exactly once";
    }
    return {a, b, c};
  }

  explicit constexpr TGeoAxes(Checked v) : m_axes{v.a, v.b, v.c, '\0'} {}

  char m_axes[4];
};

/// @}

}  // namespace ActsPlugins
