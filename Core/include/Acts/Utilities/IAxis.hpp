// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisFwd.hpp"

#include <iosfwd>
#include <vector>

namespace Acts {

/// Common base class for all Axis instance. This allows generice handling
/// such as for inspection.
class IAxis {
 public:
  /// @brief returns whether the axis is equidistant
  ///
  /// @return bool is equidistant
  virtual bool isEquidistant() const = 0;

  /// @brief returns whether the axis is variable
  ///
  /// @return bool is variable
  virtual bool isVariable() const = 0;

  /// @brief returns the type of the axis
  /// @return @c AxisType of this axis
  virtual AxisType getType() const = 0;

  /// @brief returns the boundary type set in the template param
  ///
  /// @return @c AxisBoundaryType of this axis
  virtual AxisBoundaryType getBoundaryType() const = 0;

  /// @brief Return a vector of bin edges
  /// @return Vector which contains the bin edges
  virtual std::vector<ActsScalar> getBinEdges() const = 0;

  /// @brief get minimum of binning range
  ///
  /// @return minimum of binning range
  virtual ActsScalar getMin() const = 0;

  /// @brief get maximum of binning range
  ///
  /// @return maximum of binning range
  virtual ActsScalar getMax() const = 0;

  /// @brief get total number of bins
  ///
  /// @return total number of bins (excluding under-/overflow bins)
  virtual std::size_t getNBins() const = 0;

  /// Helper function that dispatches from the @c IAxis base class
  /// to a concrete axis type. It will call the provided @p callable
  /// with a const reference to the concrete axis type.
  /// @tparam callable_t the callable type
  /// @param callable the callable object
  template <typename callable_t>
  decltype(auto) visit(const callable_t& callable) const {
    auto switchOnType = [this, &callable](auto bdtTag) {
      constexpr auto bdt = decltype(bdtTag)::value;
      switch (getType()) {
        case AxisType::Equidistant:
          return callable(
              dynamic_cast<const Axis<AxisType::Equidistant, bdt>&>(*this));
        case AxisType::Variable:
          return callable(
              dynamic_cast<const Axis<AxisType::Variable, bdt>&>(*this));
        default:
          throw std::logic_error{"Unknown type case"};
      }
    };

    switch (getBoundaryType()) {
      case AxisBoundaryType::Open:
        return switchOnType(AxisOpen);
      case AxisBoundaryType::Bound:
        return switchOnType(AxisBound);
      case AxisBoundaryType::Closed:
        return switchOnType(AxisClosed);
      default:
        throw std::logic_error{"Unknown boundary type case"};
    }
  }

  /// Check if two axes are equal
  /// @param lhs first axis
  /// @param rhs second axis
  /// @return true if the axes are equal
  friend bool operator==(const IAxis& lhs, const IAxis& rhs);

  /// Check if two axes are not equal
  /// @param lhs first axis
  /// @param rhs second axis
  /// @return true if the axes are not equal
  friend bool operator!=(const IAxis& lhs, const IAxis& rhs) {
    return !(lhs == rhs);
  }

  /// Output stream operator
  /// @param os output stream
  /// @param axis the axis to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os, const IAxis& axis);

 protected:
  /// Dispatch to the correct stream operator
  /// @param os output stream
  virtual void toStream(std::ostream& os) const = 0;
};

}  // namespace Acts
