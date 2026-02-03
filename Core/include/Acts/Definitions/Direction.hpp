// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <string>

namespace Acts {

/// The direction is always with respect to a given momentum, surface normal or
/// other general axes
class Direction final {
 private:
  enum class Value : int {
    Negative = -1,
    Positive = 1,
  };

 public:
  /// Create negative direction (-1)
  /// @return Direction with negative value
  static constexpr Direction Negative() { return Direction{Value::Negative}; }
  /// Create positive direction (+1)
  /// @return Direction with positive value
  static constexpr Direction Positive() { return Direction{Value::Positive}; }

  /// Create backward direction (equivalent to negative)
  /// @return Direction with negative value for backward propagation
  static constexpr Direction Backward() { return Direction{Value::Negative}; }
  /// Create forward direction (equivalent to positive)
  /// @return Direction with positive value for forward propagation
  static constexpr Direction Forward() { return Direction{Value::Positive}; }

  /// Create direction opposite to normal (negative)
  /// @return Direction with negative value, opposite to surface normal
  static constexpr Direction OppositeNormal() {
    return Direction{Value::Negative};
  }
  /// Create direction along normal (positive)
  /// @return Direction with positive value, along surface normal
  static constexpr Direction AlongNormal() {
    return Direction{Value::Positive};
  }

  /// This turns a signed value into a direction. Will assert on zero.
  ///
  /// @param scalar is the signed value
  ///
  /// @return a direction enum
  static constexpr Direction fromScalar(double scalar) {
    assert(scalar != 0);
    return scalar >= 0 ? Positive() : Negative();
  }

  /// This turns a signed value into a direction and 0 will be handled as a
  /// positive direction. Only use this when you are convinced that the 0 case
  /// is properly handled downstream.
  ///
  /// @param scalar is the signed value
  ///
  /// @return a direction enum
  static constexpr Direction fromScalarZeroAsPositive(double scalar) {
    return scalar >= 0 ? Positive() : Negative();
  }

  /// Convert and index [0,1] to a direction e.g. for sorting in
  /// std::array<T, 2u>
  ///
  /// @param index is the direction at input
  /// @return Direction corresponding to the index (0->Negative, 1->Positive)
  static constexpr Direction fromIndex(std::size_t index) {
    return index == 0u ? Negative() : Positive();
  }

  /// Convert dir to index [0,1] which allows to store direction dependent
  /// objects in std::array<T, 2u>
  ///
  /// @return either 0 or 1
  constexpr std::size_t index() const {
    if (m_value == Value::Negative) {
      return 0u;
    }
    return 1u;
  }

  /// Turns the direction into a signed value
  ///
  /// @return a signed value
  constexpr int sign() const { return static_cast<int>(m_value); }

  /// Reverse the direction
  ///
  /// @return an opposite direction
  constexpr Direction invert() const {
    return *this == Positive() ? Negative() : Positive();
  }

  /// Convert direction to string representation
  /// @return String representation of the direction ("positive" or "negative")
  std::string toString() const;

  /// Check if two directions are equal
  constexpr bool operator==(const Direction& rhs) const noexcept = default;

  /// Stream operator for Direction
  /// @param os Output stream
  /// @param dir Direction to output
  /// @return Reference to output stream
  friend std::ostream& operator<<(std::ostream& os, Direction dir) {
    os << dir.toString();
    return os;
  }

  /// Arithmetic operators
  /// @{

  // Direction * T

  /// Multiply Direction with integer
  /// @param dir Direction value
  /// @param value Integer to multiply
  /// @return Signed integer result
  friend constexpr int operator*(Direction dir, int value) {
    return dir.sign() * value;
  }

  /// Multiply Direction with float
  /// @param dir Direction value
  /// @param value Float to multiply
  /// @return Signed float result
  friend constexpr float operator*(Direction dir, float value) {
    return static_cast<float>(dir.sign()) * value;
  }

  /// Multiply Direction with double
  /// @param dir Direction value
  /// @param value Double to multiply
  /// @return Signed double result
  friend constexpr double operator*(Direction dir, double value) {
    return static_cast<double>(dir.sign()) * value;
  }

  /// Multiply Direction with Vector3
  /// @param dir Direction value
  /// @param value Vector3 to multiply
  /// @return Signed Vector3 result
  friend inline Acts::Vector3 operator*(Direction dir,
                                        const Acts::Vector3& value) {
    return static_cast<float>(dir.sign()) * value;
  }

  // T * Direction

  /// Multiply integer with Direction
  /// @param value Integer to multiply
  /// @param dir Direction value
  /// @return Signed integer result
  friend constexpr int operator*(int value, Direction dir) {
    return value * dir.sign();
  }

  /// Multiply float with Direction
  /// @param value Float to multiply
  /// @param dir Direction value
  /// @return Signed float result
  friend constexpr float operator*(float value, Direction dir) {
    return value * static_cast<float>(dir.sign());
  }

  /// Multiply double with Direction
  /// @param value Double to multiply
  /// @param dir Direction value
  /// @return Signed double result
  friend constexpr double operator*(double value, Direction dir) {
    return value * dir.sign();
  }

  /// Multiply Vector3 with Direction
  /// @param value Vector3 to multiply
  /// @param dir Direction value
  /// @return Signed Vector3 result
  friend Acts::Vector3 operator*(const Acts::Vector3& value, Direction dir) {
    return value * dir.sign();
  }

  // T *= Direction

  /// Multiply-assign integer with Direction
  /// @param value Integer reference to modify
  /// @param dir Direction value
  /// @return Reference to modified integer
  friend constexpr int operator*=(int& value, Direction dir) {
    value *= dir.sign();
    return value;
  }

  /// Multiply-assign float with Direction
  /// @param value Float reference to modify
  /// @param dir Direction value
  /// @return Reference to modified float
  friend constexpr float operator*=(float& value, Direction dir) {
    value *= static_cast<float>(dir.sign());
    return value;
  }

  /// Multiply-assign double with Direction
  /// @param value Double reference to modify
  /// @param dir Direction value
  /// @return Reference to modified double
  friend constexpr double operator*=(double& value, Direction dir) {
    value *= dir.sign();
    return value;
  }

  /// Multiply-assign Vector3 with Direction
  /// @param value Vector3 reference to modify
  /// @param dir Direction value
  /// @return Reference to modified Vector3
  friend Acts::Vector3& operator*=(Acts::Vector3& value, Direction dir) {
    value *= dir.sign();
    return value;
  }

  /// @}

 private:
  explicit constexpr Direction(Value value) : m_value(value) {}

  Value m_value = Value::Positive;
};

}  // namespace Acts
