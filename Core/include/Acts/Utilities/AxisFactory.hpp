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

#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace Acts {

/// @brief Per-axis resolution input supplied by the consumer of a deferred
/// axis description.
///
/// A deferred @c AxisFactory only describes the binning structure of an axis;
/// the range and the boundary type are determined by the consumer, e.g. from
/// the bounds of the surface the axis is attached to.
struct AxisResolution {
  /// Minimum edge of the axis
  double min{};
  /// Maximum edge of the axis
  double max{};
  /// Boundary type of the axis
  AxisBoundaryType boundaryType{AxisBoundaryType::Bound};
};

/// @brief Variant-like description of an axis that produces @c IAxis objects
///
/// This class captures the closed set of ways an axis can be constructed:
///
/// - @c EquidistantParams : a fully specified equidistant axis with boundary
///   type, range and number of bins.
/// - @c VariableParams : a fully specified variable axis with boundary type
///   and bin edges.
/// - @c DeferredEquidistantParams : only the number of bins is known; the
///   range and boundary type are supplied by the consumer at resolution time.
/// - @c DeferredVariableParams : only the relative bin edge distribution is
///   known, expressed as strictly increasing values normalized to [0, 1];
///   at resolution time the values are scaled affinely onto the consumer
///   supplied range.
///
/// The deferred alternatives model use cases like proto material binning,
/// where the user configures the binning structure while the axis range (and
/// whether the axis is closed) is only known once the final surface exists.
///
/// The resolution API fails fast: a fully specified axis cannot be resolved
/// with a range (there is no silent override), and a deferred description
/// cannot be resolved without one.
///
/// Optionally an @c AxisDirection can be attached, which is forwarded into
/// the created @c IAxis. A consumer that derives the direction from context
/// (e.g. from the surface) can pass its expectation to the resolution calls;
/// a mismatch with the stored direction is an error.
class AxisFactory {
 public:
  /// Parameters for a fully specified equidistant axis
  struct EquidistantParams {
    /// Boundary type of the axis
    AxisBoundaryType boundaryType{};
    /// Minimum edge of the axis
    double min{};
    /// Maximum edge of the axis
    double max{};
    /// Number of bins
    std::size_t nBins{};

    /// Check if two parameter sets are equal
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const EquidistantParams& lhs,
                           const EquidistantParams& rhs) = default;
  };

  /// Parameters for a fully specified variable axis
  struct VariableParams {
    /// Boundary type of the axis
    AxisBoundaryType boundaryType{};
    /// Bin edges, strictly increasing
    std::vector<double> edges;

    /// Check if two parameter sets are equal
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const VariableParams& lhs,
                           const VariableParams& rhs) = default;
  };

  /// Parameters for a deferred equidistant axis: only the binning structure
  /// is known, range and boundary type come from the consumer
  struct DeferredEquidistantParams {
    /// Number of bins
    std::size_t nBins{};

    /// Check if two parameter sets are equal
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const DeferredEquidistantParams& lhs,
                           const DeferredEquidistantParams& rhs) = default;
  };

  /// Parameters for a deferred variable axis: only the relative bin edge
  /// distribution is known, range and boundary type come from the consumer
  struct DeferredVariableParams {
    /// Relative bin edges, strictly increasing, with first value 0 and last
    /// value 1
    std::vector<double> normalizedEdges;

    /// Check if two parameter sets are equal
    /// @param lhs first parameter set
    /// @param rhs second parameter set
    /// @return true if the parameter sets are equal
    friend bool operator==(const DeferredVariableParams& lhs,
                           const DeferredVariableParams& rhs) = default;
  };

 private:
  /// Underlying variant type
  using Variant =
      std::variant<EquidistantParams, VariableParams, DeferredEquidistantParams,
                   DeferredVariableParams>;

  /// Construct from variant
  /// @param variant the alternative to hold
  /// @param direction the optional axis direction
  explicit AxisFactory(Variant variant, std::optional<AxisDirection> direction);

 public:
  /// Fully specified equidistant axis
  /// @param boundaryType the boundary type of the axis
  /// @param min the minimum edge of the axis
  /// @param max the maximum edge of the axis
  /// @param nBins the number of bins
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if min >= max or nBins == 0
  /// @return AxisFactory holding the equidistant description
  static AxisFactory Equidistant(
      AxisBoundaryType boundaryType, double min, double max, std::size_t nBins,
      std::optional<AxisDirection> direction = std::nullopt);

  /// Fully specified variable axis
  /// @param boundaryType the boundary type of the axis
  /// @param edges the bin edges, strictly increasing
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if fewer than two edges are given or the
  ///         edges are not strictly increasing
  /// @return AxisFactory holding the variable description
  static AxisFactory Variable(
      AxisBoundaryType boundaryType, std::vector<double> edges,
      std::optional<AxisDirection> direction = std::nullopt);

  /// Deferred equidistant axis: the range and boundary type are supplied by
  /// the consumer at resolution time
  /// @param nBins the number of bins
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if nBins == 0
  /// @return AxisFactory holding the deferred equidistant description
  static AxisFactory DeferredEquidistant(
      std::size_t nBins, std::optional<AxisDirection> direction = std::nullopt);

  /// Deferred variable axis: the normalized edges are scaled onto the range
  /// supplied by the consumer at resolution time
  /// @param normalizedEdges the relative bin edges, strictly increasing, with
  ///        first value 0 and last value 1
  /// @param direction the optional direction of the axis
  /// @throws std::invalid_argument if fewer than two values are given, the
  ///         values are not strictly increasing, or the first and last values
  ///         are not exactly 0 and 1
  /// @return AxisFactory holding the deferred variable description
  static AxisFactory DeferredVariable(
      std::vector<double> normalizedEdges,
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

  /// Get the deferred counterpart of this description: a fully specified
  /// equidistant axis keeps only its number of bins, a fully specified
  /// variable axis keeps its edges normalized to [0, 1]; deferred
  /// descriptions are returned unchanged
  /// @return AxisFactory holding the deferred description
  AxisFactory toDeferred() const;

  /// Check if the description is deferred, i.e. requires a consumer supplied
  /// range and boundary type to produce an axis
  /// @return true if the description is deferred
  bool isDeferred() const;

  /// Check if the description produces an equidistant axis
  /// @return true for the equidistant and deferred equidistant alternatives
  bool isEquidistant() const;

  /// Check if the description produces a variable axis
  /// @return true for the variable and deferred variable alternatives
  bool isVariable() const;

  /// Get the optional direction of the axis
  /// @return the direction if set
  std::optional<AxisDirection> direction() const;

  /// Get the boundary type of the axis
  /// @return the boundary type for fully specified descriptions, nullopt for
  ///         deferred descriptions
  std::optional<AxisBoundaryType> boundaryType() const;

  /// Get the number of bins
  /// @return the number of bins, defined for all alternatives
  std::size_t nBins() const;

  /// Get the description as fully specified equidistant parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the equidistant parameters
  const EquidistantParams& asEquidistant() const;

  /// Get the description as fully specified variable parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the variable parameters
  const VariableParams& asVariable() const;

  /// Get the description as deferred equidistant parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the deferred equidistant parameters
  const DeferredEquidistantParams& asDeferredEquidistant() const;

  /// Get the description as deferred variable parameters
  /// @throws std::bad_variant_access if another alternative is held
  /// @return reference to the deferred variable parameters
  const DeferredVariableParams& asDeferredVariable() const;

  /// Produce the axis from a fully specified description
  /// @param direction the direction expected by the caller; if both this and
  ///        the stored direction are set they have to agree
  /// @throws std::domain_error if the description is deferred
  /// @throws std::invalid_argument if the directions mismatch
  /// @return the created axis
  std::unique_ptr<IAxis> toAxis(
      std::optional<AxisDirection> direction = std::nullopt) const;

  /// Produce the axis from a deferred description with the consumer supplied
  /// range and boundary type
  /// @param resolution the range and boundary type to resolve with
  /// @param direction the direction expected by the caller; if both this and
  ///        the stored direction are set they have to agree
  /// @throws std::domain_error if the description is fully specified
  /// @throws std::invalid_argument if the directions mismatch or the range is
  ///         invalid
  /// @return the created axis
  std::unique_ptr<IAxis> toAxis(
      const AxisResolution& resolution,
      std::optional<AxisDirection> direction = std::nullopt) const;

  /// Get a string representation of this description
  /// @return the string representation
  std::string toString() const;

  /// Check if two descriptions are equal
  /// @param lhs first description
  /// @param rhs second description
  /// @return true if alternative, parameters and direction are equal
  friend bool operator==(const AxisFactory& lhs, const AxisFactory& rhs) {
    return lhs.m_variant == rhs.m_variant && lhs.m_direction == rhs.m_direction;
  }

  /// Output stream operator
  /// @param os output stream
  /// @param axisFactory the description to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const AxisFactory& axisFactory);

 private:
  /// Resolve the effective direction from the stored and the caller supplied
  /// direction, throwing on mismatch
  /// @param direction the caller supplied direction
  /// @return the effective direction
  std::optional<AxisDirection> resolveDirection(
      std::optional<AxisDirection> direction) const;

  Variant m_variant;
  std::optional<AxisDirection> m_direction;
};

}  // namespace Acts
