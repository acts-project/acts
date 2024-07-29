// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>

namespace Acts {

/// @enum TrackStateFlag
///
/// This enum describes the type of TrackState
enum TrackStateFlag {
  MeasurementFlag = 0,
  ParameterFlag = 1,
  OutlierFlag = 2,
  HoleFlag = 3,
  MaterialFlag = 4,
  SharedHitFlag = 5,
  NumTrackStateFlags = 6
};

class ConstTrackStateType;

/// View type over a bitset stored in a 64 bit integer
/// This view allows modifications.
class TrackStateType {
 public:
  using raw_type = std::uint64_t;
  static constexpr std::size_t kRawBits =
      std::numeric_limits<std::make_unsigned<raw_type>::type>::digits;
  /// Constructor from a reference to the underlying value container
  /// @param raw the value container
  TrackStateType(raw_type& raw) : m_raw{&raw} {}

  /// Assign the value from another set of flags
  /// @param other the other set of flags to assign
  /// @return this object
  TrackStateType& operator=(const TrackStateType& other) {
    assert(other.m_raw != nullptr);
    *m_raw = *other.m_raw;
    return *this;
  }

  /// Assign the value from another set of flags
  /// @param other the other set of flags to assign
  /// @return this object
  TrackStateType& operator=(const ConstTrackStateType& other);

  /// Automatically convert to const track state type
  operator ConstTrackStateType();

  /// Return if the bit at position @p pos is 1
  /// @param pos the bit position
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    return bs.test(pos);
  }

  /// Change the value of the bit at position @p pos to @p value.
  /// @param pos the position of the bit to change
  /// @param value the value to change the bit to
  void set(std::size_t pos, bool value = true) {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    bs.set(pos, value);
    *m_raw = bs.to_ullong();
  }

  /// Change the value of the bit at position at @p pos to @c false
  /// @param pos the position of the bit to change
  void reset(std::size_t pos) { set(pos, false); }

  friend std::ostream& operator<<(std::ostream& os, TrackStateType t) {
    assert(t.m_raw != nullptr);
    std::bitset<kRawBits> bs{*t.m_raw};
    std::bitset<TrackStateFlag::NumTrackStateFlags> trunc;
    for (std::size_t i = 0; i < TrackStateFlag::NumTrackStateFlags; i++) {
      trunc[i] = bs[i];
    }
    // SharedhitMaterialHoleOutlierParameterMeasurement
    os << "SMHOPM=" << trunc;
    return os;
  }

 private:
  raw_type* m_raw{nullptr};
};

/// View type over a bitset stored in a 64 bit integer
/// This view does not allow modifications
class ConstTrackStateType {
 public:
  using raw_type = std::uint64_t;
  static constexpr std::size_t kRawBits =
      std::numeric_limits<std::make_unsigned<raw_type>::type>::digits;

  /// Constructor from a reference to the underlying value container
  /// @param raw the value container
  ConstTrackStateType(const raw_type& raw) : m_raw{&raw} {}

  /// Return if the bit at position @p pos is 1
  /// @param pos the bit position
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const {
    assert(m_raw != nullptr);
    std::bitset<kRawBits> bs{*m_raw};
    return bs.test(pos);
  }

  friend std::ostream& operator<<(std::ostream& os, ConstTrackStateType t) {
    assert(t.m_raw != nullptr);
    std::bitset<kRawBits> bs{*t.m_raw};
    std::bitset<TrackStateFlag::NumTrackStateFlags> trunc;
    for (std::size_t i = 0; i < TrackStateFlag::NumTrackStateFlags; i++) {
      trunc[i] = bs[i];
    }
    // SharedhitMaterialHoleOutlierParameterMeasurement
    os << "SMHOPM=" << trunc;
    return os;
  }

 private:
  friend class TrackStateType;
  const raw_type* m_raw{nullptr};
};

inline TrackStateType& TrackStateType::operator=(
    const ConstTrackStateType& other) {
  assert(other.m_raw != nullptr);
  *m_raw = *other.m_raw;
  return *this;
}
inline TrackStateType::operator ConstTrackStateType() {
  return {*m_raw};
}

}  // namespace Acts
