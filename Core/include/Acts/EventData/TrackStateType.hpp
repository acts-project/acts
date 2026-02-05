// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>

namespace Acts {

/// Enumeration of track state flags.
///
/// The values define the bit positions in the underlying storage. The user may
/// define custom flags starting from bit position `NumFlags` but it is
/// recommended to start from bit position 16 or in reverse from 63 to avoid
/// conflicts with future additions.
enum class TrackStateFlag {
  MeasurementFlag [[deprecated("Replaced by HasMeasurement; consider to use "
                               "isMeasurement()/setIsMeasurement() instead")]] =
      0,
  ParameterFlag [[deprecated("Replaced by HasParameters; consider to use "
                             "hasParameters()/setHasParameters() instead")]] =
      1,
  OutlierFlag [[deprecated("Replaced by IsOutlier; consider to use "
                           "isOutlier()/setIsOutlier() instead")]] = 2,
  HoleFlag [[deprecated("Replaced by IsHole; consider to use "
                        "isHole()/setIsHole() instead")]] = 3,
  MaterialFlag [[deprecated("Replaced by HasMaterial; consider to use "
                            "hasMaterial()/setHasMaterial() instead")]] = 4,
  SharedHitFlag [[deprecated("Replaced by IsSharedHit; consider to use "
                             "isSharedHit()/setIsSharedHit() instead")]] = 5,
  SplitHitFlag [[deprecated("Replaced by IsSplitHit; consider to use "
                            "isSplitHit()/setIsSplitHit() instead")]] = 6,
  NoExpectedHitFlag
  [[deprecated("Replaced by HasNoExpectedHit; consider to use "
               "hasNoExpectedHit()/setHasNoExpectedHit() instead")]] = 7,
  NumTrackStateFlags [[deprecated("Replaced by NumFlags")]] = 8,

  /// Indicates that the track state has an associated measurement.
  /// Note that an outlier also has a measurement.
  HasMeasurement = 0,
  /// Indicates that the track state has associated parameters.
  HasParameters = 1,
  /// Indicates that the track state has an outlier measurement.
  IsOutlier = 2,
  /// Indicates that the track state has a hole (missing measurement).
  IsHole = 3,
  /// Indicates that the track state has associated material.
  HasMaterial = 4,
  /// Indicates that the track state has a shared hit measurement.
  IsSharedHit = 5,
  /// Indicates that the track state has a split hit measurement.
  IsSplitHit = 6,
  /// Indicates that the track state has no expected hit.
  HasNoExpectedHit = 7,
  /// Number of defined flags.
  NumFlags = 8
};

/// CRTP base class for @c TrackStateType and @c TrackStateTypeMap
/// @tparam Derived The derived class
/// @tparam ReadOnly Whether the derived class is read-only
template <typename Derived, bool ReadOnly>
class TrackStateTypeBase {
 public:
  /// Type alias for underlying raw data type
  using raw_type = std::uint64_t;
  /// Number of bits available in the raw storage type
  static constexpr std::size_t kRawBits =
      std::numeric_limits<std::make_unsigned_t<raw_type>>::digits;
  /// Type alias for bitset representation
  using bitset_type = std::bitset<kRawBits>;

  using enum TrackStateFlag;

  /// Assigns the flags from another TrackStateTypeBase
  /// @tparam DerivedOther The derived type of the other object
  /// @tparam ReadOnlyOther Whether the other object is read-only
  /// @param other The track state type to copy from
  /// @return Reference to this object
  template <typename DerivedOther, bool ReadOnlyOther>
  Derived operator=(
      const TrackStateTypeBase<DerivedOther, ReadOnlyOther>& other)
    requires(!ReadOnly)
  {
    self().raw() = static_cast<const DerivedOther&>(other).raw();
    return self();
  }

  /// Checks if the track state has parameters
  /// @return true if parameters are present
  bool hasParameters() const { return test(HasParameters); }

  /// Checks if the track state has material
  /// @note use isMaterial() to check for a pure material state
  /// @return true if material is present
  bool hasMaterial() const { return test(HasMaterial); }

  /// Checks if the track state **has** a measurement
  /// @note use isMeasurement() to check for a pure measurement state
  /// @return true if a measurement is present
  bool hasMeasurement() const { return test(HasMeasurement); }

  /// Checks if the track state **is** a pure material state
  /// @note use hasMaterial() to check for the presence of material
  /// @return true if it is a pure material state
  bool isMaterial() const { return test(HasMaterial) && !test(HasMeasurement); }

  /// Checks if the track state **is** a pure measurement state
  /// @note use hasMeasurement() to check for the presence of a measurement
  /// @return true if it is a pure measurement state
  bool isMeasurement() const {
    return test(HasMeasurement) && !test(IsOutlier);
  }

  /// Checks if the track state is an outlier
  /// @return true if it is an outlier
  bool isOutlier() const { return test(IsOutlier); }

  /// Checks if the track state is a hole
  /// @return true if it is a hole
  bool isHole() const { return test(IsHole); }

  /// Checks if the track state has a shared hit
  /// @return true if it has a shared hit
  bool isSharedHit() const { return test(IsSharedHit); }

  /// Checks if the track state has a split hit
  /// @return true if it has a split hit
  bool isSplitHit() const { return test(IsSplitHit); }

  /// Checks if the track state has no expected hit
  /// @return true if it has no expected hit
  bool hasNoExpectedHit() const { return test(HasNoExpectedHit); }

  /// Sets whether the track state has parameters
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setHasParameters(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasParameters, value);
    assertConsistency();
    return self();
  }

  /// Sets whether the track state has material
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setHasMaterial(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasMaterial, value);
    assertConsistency();
    return self();
  }

  /// Sets whether the track state has a measurement
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setHasMeasurement(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, value);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be a pure material state
  /// @return self-reference for chaining
  Derived& setIsMaterial()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, false);
    setUnchecked(IsOutlier, false);
    setUnchecked(IsHole, false);
    setUnchecked(IsSharedHit, false);
    setUnchecked(IsSplitHit, false);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(HasMaterial, true);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be a pure measurement state
  /// @return self-reference for chaining
  Derived& setIsMeasurement()
    requires(!ReadOnly)
  {
    setUnchecked(IsOutlier, false);
    setUnchecked(IsHole, false);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(HasMeasurement, true);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be an outlier
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setIsOutlier(bool value = true)
    requires(!ReadOnly)
  {
    if (value) {
      setUnchecked(HasMeasurement, true);
      setUnchecked(IsHole, false);
      setUnchecked(HasNoExpectedHit, false);
    }
    setUnchecked(IsOutlier, value);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be a hole
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setIsHole(bool value = true)
    requires(!ReadOnly)
  {
    if (value) {
      setUnchecked(HasMeasurement, false);
      setUnchecked(IsOutlier, false);
      setUnchecked(HasNoExpectedHit, false);
    }
    setUnchecked(IsHole, value);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be a shared hit
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setIsSharedHit(bool value = true)
    requires(!ReadOnly)
  {
    if (value) {
      setUnchecked(HasMeasurement, true);
      setUnchecked(HasNoExpectedHit, false);
    }
    setUnchecked(IsSharedHit, value);
    assertConsistency();
    return self();
  }

  /// Sets the track state to be a split hit
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setIsSplitHit(bool value = true)
    requires(!ReadOnly)
  {
    if (value) {
      setUnchecked(HasMeasurement, true);
      setUnchecked(HasNoExpectedHit, false);
    }
    setUnchecked(IsSplitHit, value);
    assertConsistency();
    return self();
  }

  /// Sets the track state to have no expected hit
  /// @param value the value to set
  /// @return self-reference for chaining
  Derived& setHasNoExpectedHit(bool value = true)
    requires(!ReadOnly)
  {
    if (value) {
      setUnchecked(HasMeasurement, false);
      setUnchecked(IsOutlier, false);
      setUnchecked(IsHole, false);
      setUnchecked(IsSharedHit, false);
      setUnchecked(IsSplitHit, false);
    }
    setUnchecked(HasNoExpectedHit, value);
    assertConsistency();
    return self();
  }

  /// Resets all flags to zero
  void reset()
    requires(!ReadOnly)
  {
    Derived::raw() = 0;
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const TrackStateTypeBase& t) {
    bool first = true;
    const auto append = [&](const char* name, bool condition) {
      if (condition) {
        if (!first) {
          os << ",";
        }
        os << name;
        first = false;
      }
    };
    append("HasParams", t.hasParameters());
    if (t.isMaterial()) {
      append("IsMat", true);
    } else {
      append("HasMat", t.hasMaterial());
    }
    if (t.isMeasurement()) {
      append("IsMeas", true);
    } else if (t.isOutlier()) {
      append("IsOutlier", true);
    } else if (t.isHole()) {
      append("IsHole", true);
    } else {
      append("HasMeas", t.hasMeasurement());
    }
    append("IsSharedHit", t.isSharedHit());
    append("IsSplitHit", t.isSplitHit());
    append("NoExpectedHit", t.hasNoExpectedHit());
    return os;
  }

  /// Return if the bit at position @p pos is set
  /// @param pos the position of the bit to test
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const { return bits().test(pos); }

  /// Return if the bit for @p flag is set
  /// @param flag the flag to test
  /// @return if the bit for @p flag is one or not
  bool test(TrackStateFlag flag) const { return test(toUnderlying(flag)); }

  /// Change the value of the bit at position @p pos to @p value.
  /// @param pos the position of the bit to change
  /// @param value the value to change the bit to
  void setUnchecked(std::size_t pos, bool value = true)
    requires(!ReadOnly)
  {
    self().raw() = bits().set(pos, value).to_ullong();
  }

  /// Change the value of the bit for @p flag to @p value.
  /// @param flag the flag to change
  /// @param value the value to change the bit to
  void setUnchecked(TrackStateFlag flag, bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(toUnderlying(flag), value);
  }

 protected:
  /// Cast to derived class
  /// @return Reference to derived object
  Derived& self() { return static_cast<Derived&>(*this); }
  /// Cast to derived class (const)
  /// @return Const reference to derived object
  const Derived& self() const { return static_cast<const Derived&>(*this); }

  /// Returns the bitset representation of the underlying raw data
  /// @return the bitset representation
  bitset_type bits() const { return bitset_type(self().raw()); }

  /// Asserts the consistency of the current flag combination
  void assertConsistency() const {
    assert(!(test(HasNoExpectedHit) &&
             (test(HasMeasurement) || test(IsOutlier) || test(IsHole) ||
              test(IsSharedHit) || test(IsSplitHit))) &&
           "TrackStateType - NoExpectedHit cannot be set with other "
           "measurement flags");
    assert(!(test(IsOutlier) && test(IsHole)) &&
           "TrackStateType - Outlier and Hole cannot be set simultaneously");
    assert(
        !(test(IsHole) && test(HasMeasurement)) &&
        "TrackStateType - Hole and Measurement cannot be set simultaneously");
    assert(!(test(IsOutlier) && !test(HasMeasurement)) &&
           "TrackStateType - Outlier flag requires Measurement flag to be set");
  }
};

/// @c TrackStateTypeBase captured by value.
class TrackStateType : public TrackStateTypeBase<TrackStateType, false> {
 public:
  /// Base class type
  using Base = TrackStateTypeBase<TrackStateType, false>;

  TrackStateType() = default;
  /// Construct from raw value
  /// @param raw The raw flag bits
  explicit TrackStateType(raw_type raw) : m_raw{raw} {
    Base::assertConsistency();
  }

  using Base::operator=;

  /// Access raw storage
  /// @return Reference to raw storage
  raw_type& raw() { return m_raw; }
  /// Access raw storage (const)
  /// @return Const reference to raw storage
  const raw_type& raw() const { return m_raw; }

 private:
  raw_type m_raw{0};
};

/// @c TrackStateTypeBase mapped to external storage.
/// @tparam ReadOnly Whether the mapped storage is read-only
template <bool ReadOnly>
class TrackStateTypeMap
    : public TrackStateTypeBase<TrackStateTypeMap<ReadOnly>, ReadOnly> {
 public:
  /// Base class type
  using Base = TrackStateTypeBase<TrackStateTypeMap<ReadOnly>, ReadOnly>;
  /// Underlying raw storage type, const-qualified based on ReadOnly
  using raw_type = const_if_t<ReadOnly, typename Base::raw_type>;

  /// Constructor from external raw storage reference
  /// @param raw_ref Reference to external raw storage
  explicit TrackStateTypeMap(raw_type& raw_ref) : m_raw_ptr{&raw_ref} {
    assert(m_raw_ptr != nullptr && "TrackStateTypeMap - raw reference is null");
    Base::assertConsistency();
  }

  using Base::operator=;

  /// Access raw storage
  /// @return Reference to raw storage
  raw_type& raw()
    requires(!ReadOnly)
  {
    return *m_raw_ptr;
  }
  /// Access raw storage (const)
  /// @return Const reference to raw storage
  const raw_type& raw() const { return *m_raw_ptr; }

 private:
  raw_type* m_raw_ptr{nullptr};
};

/// Mutable track state type map allowing modification
using MutableTrackStateTypeMap = TrackStateTypeMap<false>;
/// Const track state type map for read-only access
using ConstTrackStateTypeMap = TrackStateTypeMap<true>;

}  // namespace Acts
