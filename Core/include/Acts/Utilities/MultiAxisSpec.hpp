// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"

#include <array>
#include <memory>
#include <ostream>
#include <span>
#include <string>
#include <vector>

namespace Acts {

class Surface;

/// @brief Specification of a multi-dimensional binning that builds the axes of
/// a grid
///
/// Bundles one @c AxisSpec per grid dimension, mirroring the
/// @c IMultiAxis / @c IMultiAxisXD split: the dimension is known at runtime
/// here and at compile time in the derived @c MultiAxisSpecXD . Unlike
/// @c IMultiAxis this is a value type and can be held by configuration
/// objects.
///
/// The build API mirrors @c AxisSpec : fully specified specs build
/// without input, deferred ones need one @c AxisSpec::Options per axis, in
/// storage order.
class MultiAxisSpec {
 public:
  /// Consumer supplied inputs, one @c AxisSpec::Options per axis, in axis
  /// order
  using Options = std::vector<AxisSpec::Options>;

  /// Construct from one axis spec per dimension
  /// @param axisSpecs the axis specs, in axis order
  /// @throws std::invalid_argument if no axis spec is given
  explicit MultiAxisSpec(std::vector<AxisSpec> axisSpecs);

  /// Copy constructor
  /// @param other the spec to copy from
  MultiAxisSpec(const MultiAxisSpec& other) = default;
  /// Move constructor
  /// @param other the spec to move from
  MultiAxisSpec(MultiAxisSpec&& other) noexcept = default;
  /// Copy assignment
  /// @param other the spec to copy from
  /// @return reference to this
  MultiAxisSpec& operator=(const MultiAxisSpec& other) = default;
  /// Move assignment
  /// @param other the spec to move from
  /// @return reference to this
  MultiAxisSpec& operator=(MultiAxisSpec&& other) noexcept = default;

  virtual ~MultiAxisSpec() = default;

  /// Get the number of axes spanning the grid
  /// @return number of axes (i.e. the dimension of the grid)
  std::size_t size() const;

  /// Get the axis spec at the given dimension
  /// @param i index of the axis
  /// @return const reference to the requested axis spec
  const AxisSpec& axisSpec(std::size_t i) const;

  /// Get all axis specs
  /// @return view onto the axis specs, in axis order
  std::span<const AxisSpec> axisSpecs() const;

  /// Check if any of the contained specs is deferred, i.e. requires
  /// consumer supplied options to produce axes
  /// @return true if any axis spec is deferred
  bool isDeferred() const;

  /// Build a multi-axis, one option set per axis or none at all
  /// @param options one option set per axis, in axis order, or empty
  /// @throws std::domain_error if a property is given by neither side, or the
  ///         dimension exceeds the supported maximum of 3
  /// @throws std::invalid_argument if the number of option sets is neither
  ///         zero nor the number of axes, or a property mismatches
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> buildMultiAxis(const Options& options = {}) const;

 protected:
  /// Build a multi-axis from a view onto the option sets
  /// @param options one option set per axis, in axis order, or empty
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> buildMultiAxisImpl(
      std::span<const AxisSpec::Options> options) const;

 public:
  /// Get a string representation of this spec
  /// @return the string representation
  std::string toString() const;

  /// Check if two specs are equal
  /// @param lhs first spec
  /// @param rhs second spec
  /// @return true if all axis specs are equal
  friend bool operator==(const MultiAxisSpec& lhs,
                         const MultiAxisSpec& rhs) = default;

  /// Output stream operator
  /// @param os output stream
  /// @param multiAxisSpec the spec to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const MultiAxisSpec& multiAxisSpec) {
    return os << multiAxisSpec.toString();
  }

 private:
  std::vector<AxisSpec> m_axisSpecs;
};

/// @brief Multi-dimensional binning spec of a compile-time dimension
///
/// Adds a statically sized construction and build API on top of
/// @c MultiAxisSpec , mirroring @c IMultiAxis and @c IMultiAxisXD .
///
/// @tparam DIM number of axes (dimension of the grid)
template <std::size_t DIM>
class MultiAxisSpecXD : public MultiAxisSpec {
 public:
  static_assert(DIM >= 1 && DIM <= 3,
                "MultiAxisSpecXD supports 1 to 3 dimensions");

  /// Construct from one axis spec per dimension
  /// @param axisSpecs the axis specs, in axis order
  explicit MultiAxisSpecXD(std::array<AxisSpec, DIM> axisSpecs)
      : MultiAxisSpec(
            std::vector<AxisSpec>(std::make_move_iterator(axisSpecs.begin()),
                                  std::make_move_iterator(axisSpecs.end()))) {}

  /// Consumer supplied inputs, one @c AxisSpec::Options per axis, in axis
  /// order; an all-default set leaves every property to the spec
  using Options = std::array<AxisSpec::Options, DIM>;

  /// Build a multi-axis, one option set per axis
  /// @param options one option set per axis, in axis order
  /// @throws std::domain_error if a property is given by neither side
  /// @throws std::invalid_argument if a property mismatches
  /// @return the created multi-axis of dimension @c DIM
  std::unique_ptr<IMultiAxisXD<DIM>> buildMultiAxis(
      const Options& options = {}) const {
    return downcast(buildMultiAxisImpl(options));
  }

 private:
  /// Downcast a runtime-dimension multi-axis to the compile-time dimension
  /// @param multiAxis the multi-axis to downcast
  /// @return the downcasted multi-axis
  static std::unique_ptr<IMultiAxisXD<DIM>> downcast(
      std::unique_ptr<IMultiAxis> multiAxis) {
    auto* xd = dynamic_cast<IMultiAxisXD<DIM>*>(multiAxis.get());
    if (xd == nullptr) {
      throw std::logic_error(
          "MultiAxisSpecXD: unexpected multi-axis dimension");
    }
    multiAxis.release();
    return std::unique_ptr<IMultiAxisXD<DIM>>(xd);
  }
};

/// Type alias for a multi-axis spec of dimension 1
using MultiAxisSpec1D = MultiAxisSpecXD<1>;
/// Type alias for a multi-axis spec of dimension 2
using MultiAxisSpec2D = MultiAxisSpecXD<2>;

/// @brief Build the multi-axis of a surface binning against that surface
///
/// A surface spans two local coordinates, so the binning is strictly
/// two-dimensional; use a single bin in one direction to bin along the other
/// only. Both specs are matched to the canonical local axis directions
/// of the surface (see @c Surface::localAxes ), either positionally if
/// neither carries a direction, or by direction if both do. Mixing the two is
/// rejected. The axes are returned in canonical direction order.
///
/// @param multiAxisSpec the binning spec to build from
/// @param surface the surface to resolve the ranges against
/// @throws std::invalid_argument if the directions cannot be matched or the
///         surface is unsupported
/// @throws std::domain_error if any spec is fully specified instead
///         of deferred
/// @return the created multi-axis, in canonical direction order
std::unique_ptr<IMultiAxis2D> resolveMultiAxis(
    const MultiAxisSpec2D& multiAxisSpec, const Surface& surface);

}  // namespace Acts
