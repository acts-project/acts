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

enum class TrackStateFlag {
  MeasurementFlag
  [[deprecated("Replaced by HasMeasurement; but you may want to use "
               "isMeasurement()/setIsMeasurement() instead")]] = 0,
  ParameterFlag
  [[deprecated("Replaced by HasParameters; but you may want to use "
               "hasParameters()/setHasParameters() instead")]] = 1,
  OutlierFlag [[deprecated("Replaced by IsOutlier; but you may want to use "
                           "isOutlier()/setIsOutlier() instead")]] = 2,
  HoleFlag [[deprecated("Replaced by IsHole; but you may want to use "
                        "isHole()/setIsHole() instead")]] = 3,
  MaterialFlag [[deprecated("Replaced by HasMaterial; but you may want to use "
                            "hasMaterial()/setHasMaterial() instead")]] = 4,
  SharedHitFlag [[deprecated("Replaced by IsSharedHit; but you may want to use "
                             "isSharedHit()/setIsSharedHit() instead")]] = 5,
  SplitHitFlag [[deprecated("Replaced by IsSplitHit; but you may want to use "
                            "isSplitHit()/setIsSplitHit() instead")]] = 6,
  NoExpectedHitFlag
  [[deprecated("Replaced by HasNoExpectedHit; but you may want to use "
               "hasNoExpectedHit()/setHasNoExpectedHit() instead")]] = 7,
  NumTrackStateFlags [[deprecated("Replaced by NumFlags")]] = 8,

  HasMeasurement = 0,
  HasParameters = 1,
  IsOutlier = 2,
  IsHole = 3,
  HasMaterial = 4,
  IsSharedHit = 5,
  IsSplitHit = 6,
  HasNoExpectedHit = 7,
  NumFlags = 8
};

template <typename Derived, bool ReadOnly = false>
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

  template <typename DerivedOther, bool ReadOnlyOther>
  Derived operator=(
      const TrackStateTypeBase<DerivedOther, ReadOnlyOther>& other)
    requires(!ReadOnly)
  {
    self().raw() = static_cast<const DerivedOther&>(other).raw();
    return self();
  }

  bool hasParameters() const { return test(HasParameters); }
  bool hasMaterial() const { return test(HasMaterial); }
  bool hasMeasurement() const { return test(HasMeasurement); }
  bool isMaterial() const { return test(HasMaterial) && !test(HasMeasurement); }
  bool isMeasurement() const {
    return test(HasMeasurement) && !test(IsOutlier);
  }
  bool isOutlier() const { return test(IsOutlier); }
  bool isHole() const { return test(IsHole); }
  bool isSharedHit() const { return test(IsSharedHit); }
  bool isSplitHit() const { return test(IsSplitHit); }
  bool hasNoExpectedHit() const { return test(HasNoExpectedHit); }
  Derived& setHasParameters(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasParameters, value);
    assertConsistency();
    return self();
  }
  Derived& setHasMaterial(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasMaterial, value);
    assertConsistency();
    return self();
  }
  Derived& setHasMeasurement(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, value);
    assertConsistency();
    return self();
  }
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
  Derived& setIsOutlier()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, true);
    setUnchecked(IsHole, false);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(IsOutlier, true);
    assertConsistency();
    return self();
  }
  Derived& setIsHole()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, false);
    setUnchecked(IsOutlier, false);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(IsHole, true);
    assertConsistency();
    return self();
  }
  Derived& setIsSharedHit()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, true);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(IsSharedHit, true);
    assertConsistency();
    return self();
  }
  Derived& setIsSplitHit()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, true);
    setUnchecked(HasNoExpectedHit, false);
    setUnchecked(IsSplitHit, true);
    assertConsistency();
    return self();
  }
  Derived& setHasNoExpectedHit()
    requires(!ReadOnly)
  {
    setUnchecked(HasMeasurement, false);
    setUnchecked(IsOutlier, false);
    setUnchecked(IsHole, false);
    setUnchecked(IsSharedHit, false);
    setUnchecked(IsSplitHit, false);
    setUnchecked(HasNoExpectedHit, true);
    assertConsistency();
    return self();
  }

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

  /// Return if the bit at position @p pos is 1
  /// @param pos the position of the bit to test
  /// @return if the bit at @p pos is one or not
  bool test(std::size_t pos) const { return bits().test(pos); }

  bool test(TrackStateFlag flag) const { return test(toUnderlying(flag)); }

  /// Change the value of the bit at position @p pos to @p value.
  /// @param pos the position of the bit to change
  /// @param value the value to change the bit to
  void setUnchecked(std::size_t pos, bool value = true)
    requires(!ReadOnly)
  {
    self().raw() = bits().set(pos, value).to_ullong();
  }

  void setUnchecked(TrackStateFlag flag, bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(toUnderlying(flag), value);
  }

 protected:
  Derived& self() { return static_cast<Derived&>(*this); }
  const Derived& self() const { return static_cast<const Derived&>(*this); }

  bitset_type bits() const { return bitset_type(self().raw()); }

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

class TrackStateType : public TrackStateTypeBase<TrackStateType> {
 public:
  using Base = TrackStateTypeBase<TrackStateType>;

  TrackStateType() = default;
  explicit TrackStateType(raw_type raw) : m_raw{raw} {
    Base::assertConsistency();
  }

  using Base::operator=;

  raw_type& raw() { return m_raw; }
  const raw_type& raw() const { return m_raw; }

 private:
  raw_type m_raw{0};
};

template <bool ReadOnly>
class TrackStateTypeMap
    : public TrackStateTypeBase<TrackStateTypeMap<ReadOnly>, ReadOnly> {
 public:
  using Base = TrackStateTypeBase<TrackStateTypeMap<ReadOnly>, ReadOnly>;
  using raw_type = const_if_t<ReadOnly, typename Base::raw_type>;

  explicit TrackStateTypeMap(raw_type& raw_ref) : m_raw_ptr{&raw_ref} {
    assert(m_raw_ptr != nullptr && "TrackStateTypeMap - raw reference is null");
    Base::assertConsistency();
  }

  using Base::operator=;

  raw_type& raw()
    requires(!ReadOnly)
  {
    return *m_raw_ptr;
  }
  const raw_type& raw() const { return *m_raw_ptr; }

 private:
  raw_type* m_raw_ptr{nullptr};
};

using MutableTrackStateTypeMap = TrackStateTypeMap<false>;
using ConstTrackStateTypeMap = TrackStateTypeMap<true>;

}  // namespace Acts
