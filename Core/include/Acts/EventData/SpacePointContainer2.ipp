// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"

#include "Acts/EventData/SpacePointProxy2.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts::Experimental {

inline MutableSpacePointProxy2 SpacePointContainer2::at(Index index) {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in SpacePointContainer2: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return MutableProxy(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::at(Index index) const {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in SpacePointContainer2: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return ConstProxy(*this, index);
}

inline MutableSpacePointProxy2 SpacePointContainer2::operator[](
    Index index) noexcept {
  return MutableProxy(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::operator[](
    Index index) const noexcept {
  return ConstProxy(*this, index);
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantX(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::X)) {
    return this->x(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XY)) {
    return this->xy(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[0];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantX: Invalid column for X coordinate");
  }
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantY(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::Y)) {
    return this->y(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XY)) {
    return this->xy(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[1];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantY: Invalid column for Y coordinate");
  }
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantZ(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::Z)) {
    return this->z(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::ZR)) {
    return this->zr(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[2];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[2];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantZ: Invalid column for Z coordinate");
  }
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantR(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::R)) {
    return this->r(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::ZR)) {
    return this->zr(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[3];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantR: Invalid column for R coordinate");
  }
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantVarianceZ(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZ)) {
    return this->varianceZ(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZR)) {
    return this->varianceZR(index)[0];
  } else {
    static_assert(false,
                  "SpacePointContainer2::variantVarianceZ: Invalid column for "
                  "VarianceZ coordinate");
  }
}

template <SpacePointColumns column>
inline float &SpacePointContainer2::variantVarianceR(Index index) noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceR)) {
    return this->varianceR(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZR)) {
    return this->varianceZR(index)[1];
  } else {
    static_assert(false,
                  "SpacePointContainer2::variantVarianceR: Invalid column for "
                  "VarianceR coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantX(Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::X)) {
    return this->x(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XY)) {
    return this->xy(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[0];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantX: Invalid column for X coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantY(Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::Y)) {
    return this->y(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XY)) {
    return this->xy(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[1];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantY: Invalid column for Y coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantZ(Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::Z)) {
    return this->z(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::ZR)) {
    return this->zr(index)[0];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZ)) {
    return this->xyz(index)[2];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[2];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantZ: Invalid column for Z coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantR(Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::R)) {
    return this->r(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::ZR)) {
    return this->zr(index)[1];
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::XYZR)) {
    return this->xyzr(index)[3];
  } else {
    static_assert(
        false,
        "SpacePointContainer2::variantR: Invalid column for R coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantVarianceZ(
    Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZ)) {
    return this->varianceZ(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZR)) {
    return this->varianceZR(index)[0];
  } else {
    static_assert(false,
                  "SpacePointContainer2::variantVarianceZ: Invalid column for "
                  "VarianceZ coordinate");
  }
}

template <SpacePointColumns column>
inline float SpacePointContainer2::variantVarianceR(
    Index index) const noexcept {
  if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceR)) {
    return this->varianceR(index);
  } else if constexpr (ACTS_CHECK_BIT(column, SpacePointColumns::VarianceZR)) {
    return this->varianceZR(index)[1];
  } else {
    static_assert(false,
                  "SpacePointContainer2::variantVarianceR: Invalid column for "
                  "VarianceR coordinate");
  }
}

}  // namespace Acts::Experimental
