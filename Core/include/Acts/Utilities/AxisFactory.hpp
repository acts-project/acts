// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/IAxis.hpp"

#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <variant>
#include <vector>

namespace Acts {

/// @brief Axis properties supplied by the consumer of an @c AxisFactory at
/// build time
///
/// Each property fills in the corresponding one of the description if that is
/// unset, and validates it otherwise. Spelled @c AxisFactory::Options at the
/// use sites; it only lives at namespace scope so that it is complete where
/// @c AxisFactory::buildAxis defaults it.
struct AxisFactoryOptions {
  /// Minimum edge of the axis
  std::optional<double> min = std::nullopt;
  /// Maximum edge of the axis
  std::optional<double> max = std::nullopt;
  /// Boundary type of the axis
  std::optional<AxisBoundaryType> boundaryType = std::nullopt;
  /// Direction of the axis
  std::optional<AxisDirection> direction = std::nullopt;
};

/// @brief Variant-like description of an axis that builds @c IAxis objects
///
/// Every property except the binning structure itself is optional. What the
/// description leaves open is supplied by the consumer as @c Options at build
/// time, e.g. from the bounds of the surface the axis is attached to. This
/// models proto material binning, where a configuration may fix the number of
/// bins but not the range, or a range but not whether the axis wraps.
///
/// Per property the rule is the same: it has to be given by exactly one side,
/// or by both with the same value. Supplying a property that the description
/// already fixes therefore validates it instead of overriding it, and a
/// property that neither side gives is an error.
class AxisFactory {
 public:
  /// Axis properties supplied by the consumer at build time
  using Options = AxisFactoryOptions;

  /// Parameters for an equidistant axis
  struct EquidistantParams {
    /// Number of bins
    std::size_t nBins = 0;
    /// Minimum edge of the axis
    std::optional<double> min = std::nullopt;
    /// Maximum edge of the axis
    std::optional<double> max = std::nullopt;

    /// Check if two parameter sets are equal, required for comparing the
    /// enclosing description
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const EquidistantParams& lhs,
                           const EquidistantParams& rhs) = default;
  };

  /// Parameters for a variable axis with absolute bin edges
  struct VariableParams {
    /// Bin edges, strictly increasing
    std::vector<double> edges;

    /// Check if two parameter sets are equal, required for comparing the
    /// enclosing description
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const VariableParams& lhs,
                           const VariableParams& rhs) = default;
  };

  /// Parameters for a variable axis whose bin edges are relative and scaled
  /// onto the range supplied by the consumer
  struct DeferredVariableParams {
    /// Relative bin edges, strictly increasing, with first value 0 and last
    /// value 1
    std::vector<double> normalizedEdges;

    /// Check if two parameter sets are equal, required for comparing the
    /// enclosing description
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const DeferredVariableParams& lhs,
                           const DeferredVariableParams& rhs) = default;
  };

 private:
  /// Underlying variant type
  using Variant =
      std::variant<EquidistantParams, VariableParams, DeferredVariableParams>;

  /// Construct from variant
  /// @param variant the alternative to hold
  /// @param boundaryType the optional axis boundary type
  /// @param direction the optional axis direction
  explicit AxisFactory(Variant variant,
                       std::optional<AxisBoundaryType> boundaryType,
                       std::optional<AxisDirection> direction);

 public:
  /// Equidistant axis; every property but the number of bins may be left to
  /// the consumer
  /// @param nBins the number of bins
  /// @param min the optional minimum edge of the axis
  /// @param max the optional maximum edge of the axis
  /// @param boundaryType the optional boundary type of the axis
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if min >= max or nBins == 0
  /// @return AxisFactory holding the equidistant description
  static AxisFactory Equidistant(
      std::size_t nBins, std::optional<double> min = std::nullopt,
      std::optional<double> max = std::nullopt,
      std::optional<AxisBoundaryType> boundaryType = std::nullopt,
      std::optional<AxisDirection> direction = std::nullopt);

  /// Equidistant axis with only the number of bins fixed
  /// @param nBins the number of bins
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if nBins == 0
  /// @return AxisFactory holding the equidistant description
  static AxisFactory DeferredEquidistant(
      std::size_t nBins, std::optional<AxisDirection> direction = std::nullopt);

  /// Variable axis with absolute bin edges
  /// @param edges the bin edges, strictly increasing
  /// @param boundaryType the optional boundary type of the axis
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if fewer than two edges are given or the
  ///         edges are not strictly increasing
  /// @return AxisFactory holding the variable description
  static AxisFactory Variable(
      std::vector<double> edges,
      std::optional<AxisBoundaryType> boundaryType = std::nullopt,
      std::optional<AxisDirection> direction = std::nullopt);

  /// Variable axis whose normalized edges are scaled onto the range supplied
  /// by the consumer
  /// @param normalizedEdges the relative bin edges, strictly increasing, with
  ///        first value 0 and last value 1
  /// @param boundaryType the optional boundary type of the axis
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if fewer than two values are given, the
  ///         values are not strictly increasing, or the first and last values
  ///         are not exactly 0 and 1
  /// @return AxisFactory holding the deferred variable description
  static AxisFactory DeferredVariable(
      std::vector<double> normalizedEdges,
      std::optional<AxisBoundaryType> boundaryType = std::nullopt,
      std::optional<AxisDirection> direction = std::nullopt);

  /// Capture an existing axis as a fully specified description
  /// @param axis the axis to decompose
  /// @return AxisFactory holding the equidistant or variable description of
  ///         the given axis, including its direction if set
  static AxisFactory FromAxis(const IAxis& axis);

  /// Get a copy of this description with the given direction attached
  /// @param direction the direction to attach
  /// @return AxisFactory with the direction set
  AxisFactory withDirection(AxisDirection direction) const;

  /// Get the counterpart that leaves every property to the consumer: an
  /// equidistant description keeps only its number of bins, a variable one
  /// its edges normalized to [0, 1]
  /// @return AxisFactory holding the deferred description
  AxisFactory toDeferred() const;

  /// Check if the description needs @c Options to build, i.e. leaves at least
  /// one property to the consumer
  /// @return true if the description is deferred
  bool isDeferred() const;

  /// Check if the description produces an equidistant axis
  /// @return true for the equidistant alternative
  bool isEquidistant() const;

  /// Check if the description produces a variable axis
  /// @return true for the two variable alternatives
  bool isVariable() const;

  /// Check if the description holds normalized instead of absolute bin edges
  /// @return true for the deferred variable alternative
  bool isDeferredVariable() const;

  /// Get the number of bins
  /// @return the number of bins, defined for all alternatives
  std::size_t nBins() const;

  /// Get the boundary type of the axis
  /// @return the boundary type if the description fixes it
  std::optional<AxisBoundaryType> boundaryType() const;

  /// Get the direction of the axis
  /// @return the direction if the description fixes it
  std::optional<AxisDirection> direction() const;

  /// Get the description as equidistant parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the equidistant parameters
  const EquidistantParams& asEquidistant() const;

  /// Get the description as variable parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the variable parameters
  const VariableParams& asVariable() const;

  /// Get the description as deferred variable parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the deferred variable parameters
  const DeferredVariableParams& asDeferredVariable() const;

  /// Build the axis, filling in the properties the description leaves open
  /// and validating the ones it fixes
  /// @param options the properties supplied by the consumer
  /// @throws std::domain_error if a property is given by neither side
  /// @throws std::invalid_argument if a property is given by both sides with
  ///         different values, or the resulting range is invalid
  /// @return the created axis
  std::unique_ptr<IAxis> buildAxis(const Options& options = {}) const;

  /// Get a string representation of this description
  /// @return the string representation
  std::string toString() const;

  /// Check if two descriptions are equal
  /// @param lhs first description
  /// @param rhs second description
  /// @return true if alternative, parameters and direction are equal
  friend bool operator==(const AxisFactory& lhs,
                         const AxisFactory& rhs) = default;

  /// Output stream operator
  /// @param os output stream
  /// @param axisFactory the description to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const AxisFactory& axisFactory) {
    return os << axisFactory.toString();
  }

 private:
  Variant m_variant;
  std::optional<AxisBoundaryType> m_boundaryType;
  std::optional<AxisDirection> m_direction;
};

}  // namespace Acts
