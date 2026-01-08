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
  Parameter = 0,
  Material = 1,
  Measurement = 2,
  Outlier = 3,
  Hole = 4,
  SharedHit = 5,
  SplitHit = 6,
  NoExpectedHit = 7,
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

  bool hasParameters() const { return test(Parameter); }
  bool hasMaterial() const { return test(Material); }
  bool hasMeasurement() const { return test(Measurement); }
  bool isMaterial() const { return test(Material) && !test(Measurement); }
  bool isMeasurement() const { return test(Measurement) && !test(Outlier); }
  bool isOutlier() const { return test(Outlier); }
  bool isHole() const { return test(Hole); }
  bool isSharedHit() const { return test(SharedHit); }
  bool isSplitHit() const { return test(SplitHit); }
  bool hasNoExpectedHit() const { return test(NoExpectedHit); }

  Derived& setHasParameters(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(Parameter, value);
    validate();
    return self();
  }
  Derived& setHasMaterial(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(Material, value);
    validate();
    return self();
  }
  Derived& setHasMeasurement(bool value = true)
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, value);
    validate();
    return self();
  }
  Derived& setIsMaterial()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, false);
    setUnchecked(Outlier, false);
    setUnchecked(Hole, false);
    setUnchecked(SharedHit, false);
    setUnchecked(SplitHit, false);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(Material, true);
    validate();
    return self();
  }
  Derived& setIsMeasurement()
    requires(!ReadOnly)
  {
    setUnchecked(Outlier, false);
    setUnchecked(Hole, false);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(Measurement, true);
    validate();
    return self();
  }
  Derived& setIsOutlier()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, true);
    setUnchecked(Hole, false);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(Outlier, true);
    validate();
    return self();
  }
  Derived& setIsHole()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, false);
    setUnchecked(Outlier, false);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(Hole, true);
    validate();
    return self();
  }
  Derived& setIsSharedHit()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, true);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(SharedHit, true);
    validate();
    return self();
  }
  Derived& setIsSplitHit()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, true);
    setUnchecked(NoExpectedHit, false);
    setUnchecked(SplitHit, true);
    validate();
    return self();
  }
  Derived& setHasNoExpectedHit()
    requires(!ReadOnly)
  {
    setUnchecked(Measurement, false);
    setUnchecked(Outlier, false);
    setUnchecked(Hole, false);
    setUnchecked(SharedHit, false);
    setUnchecked(SplitHit, false);
    setUnchecked(NoExpectedHit, true);
    validate();
    return self();
  }

  void reset()
    requires(!ReadOnly)
  {
    Derived::raw() = 0;
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const TrackStateTypeBase& t) {
    os << "TrackStateType[";
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
    append("HasParameters", t.hasParameters());
    if (t.isMaterial()) {
      append("IsMaterial", true);
    } else {
      append("HasMaterial", t.hasMaterial());
    }
    if (t.isMeasurement()) {
      append("IsMeasurement", true);
    } else if (t.isOutlier()) {
      append("IsOutlier", true);
    } else if (t.isHole()) {
      append("IsHole", true);
    } else {
      append("HasMeasurement", t.hasMeasurement());
    }
    append("IsSharedHit", t.isSharedHit());
    append("IsSplitHit", t.isSplitHit());
    append("NoExpectedHit", t.hasNoExpectedHit());
    os << "]";
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

  void validate() const {
    assert(!(test(NoExpectedHit) &&
             (test(Measurement) || test(Outlier) || test(Hole) ||
              test(SharedHit) || test(SplitHit))) &&
           "TrackStateType - NoExpectedHit cannot be set with other "
           "measurement flags");
    assert(!(test(Outlier) && test(Hole)) &&
           "TrackStateType - Outlier and Hole cannot be set simultaneously");
    assert(
        !(test(Hole) && test(Measurement)) &&
        "TrackStateType - Hole and Measurement cannot be set simultaneously");
    assert(!(test(Outlier) && !test(Measurement)) &&
           "TrackStateType - Outlier flag requires Measurement flag to be set");
  }
};

class TrackStateType : public TrackStateTypeBase<TrackStateType> {
 public:
  using Base = TrackStateTypeBase<TrackStateType>;

  TrackStateType() = default;
  explicit TrackStateType(raw_type raw) : m_raw{raw} { Base::validate(); }

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
    Base::validate();
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
