// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/IMultiAxis.hpp"

#include <array>
#include <memory>
#include <ostream>
#include <span>
#include <string>
#include <vector>

namespace Acts {

class Surface;

/// @brief Description of a multi-dimensional binning that builds the axes of
/// a grid
///
/// Bundles one @c AxisFactory per grid dimension, mirroring the
/// @c IMultiAxis / @c IMultiAxisXD split: the dimension is known at runtime
/// here and at compile time in the derived @c MultiAxisFactoryXD . Unlike
/// @c IMultiAxis this is a value type and can be held by configuration
/// objects.
///
/// The build API mirrors @c AxisFactory : fully specified descriptions build
/// without input, deferred ones need one @c AxisFactory::Options per axis, in
/// storage order.
class MultiAxisFactory {
 public:
  /// Non-owning view onto the consumer supplied inputs, one
  /// @c AxisFactory::Options per axis, in axis order
  using OptionsView = std::span<const AxisFactory::Options>;

  /// Construct from one axis description per dimension
  /// @param axisFactories the axis descriptions, in axis order
  /// @throws std::invalid_argument if no axis description is given
  explicit MultiAxisFactory(std::vector<AxisFactory> axisFactories);

  /// Copy constructor
  /// @param other the description to copy from
  MultiAxisFactory(const MultiAxisFactory& other) = default;
  /// Move constructor
  /// @param other the description to move from
  MultiAxisFactory(MultiAxisFactory&& other) noexcept = default;
  /// Copy assignment
  /// @param other the description to copy from
  /// @return reference to this
  MultiAxisFactory& operator=(const MultiAxisFactory& other) = default;
  /// Move assignment
  /// @param other the description to move from
  /// @return reference to this
  MultiAxisFactory& operator=(MultiAxisFactory&& other) noexcept = default;

  virtual ~MultiAxisFactory() = default;

  /// Get the number of axes spanning the grid
  /// @return number of axes (i.e. the dimension of the grid)
  std::size_t size() const;

  /// Get the axis description at the given dimension
  /// @param i index of the axis
  /// @return const reference to the requested axis description
  const AxisFactory& axisFactory(std::size_t i) const;

  /// Get all axis descriptions
  /// @return const reference to the axis descriptions, in axis order
  const std::vector<AxisFactory>& axisFactories() const;

  /// Check if any of the contained descriptions is deferred, i.e. requires
  /// consumer supplied options to produce axes
  /// @return true if any axis description is deferred
  bool isDeferred() const;

  /// Build the axes from fully specified descriptions
  /// @throws std::domain_error if any description is deferred
  /// @return the created axes, in axis order
  std::vector<std::unique_ptr<IAxis>> buildAxes() const;

  /// Build the axes from deferred descriptions with the consumer supplied
  /// ranges, boundary types and directions
  /// @param options one option set per axis, in axis order
  /// @throws std::domain_error if any description is fully specified
  /// @throws std::invalid_argument if the number of option sets does not
  ///         match the number of axes, or a direction mismatches
  /// @return the created axes, in axis order
  std::vector<std::unique_ptr<IAxis>> buildAxes(OptionsView options) const;

  /// Build a multi-axis from fully specified descriptions
  /// @throws std::domain_error if any description is deferred or the
  ///         dimension exceeds the supported maximum of 3
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> buildMultiAxis() const;

  /// Build a multi-axis from deferred descriptions with the consumer supplied
  /// ranges, boundary types and directions
  /// @param options one option set per axis, in axis order
  /// @throws std::domain_error if any description is fully specified or the
  ///         dimension exceeds the supported maximum of 3
  /// @throws std::invalid_argument if the number of option sets does not
  ///         match the number of axes, or a direction mismatches
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> buildMultiAxis(OptionsView options) const;

  /// Get a string representation of this description
  /// @return the string representation
  std::string toString() const;

  /// Check if two descriptions are equal
  /// @param lhs first description
  /// @param rhs second description
  /// @return true if all axis descriptions are equal
  friend bool operator==(const MultiAxisFactory& lhs,
                         const MultiAxisFactory& rhs) = default;

  /// Output stream operator
  /// @param os output stream
  /// @param multiAxisFactory the description to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const MultiAxisFactory& multiAxisFactory) {
    return os << multiAxisFactory.toString();
  }

 private:
  /// Assemble a multi-axis from already built axes
  /// @param axes the axes, in axis order
  /// @return the created multi-axis
  static std::unique_ptr<IMultiAxis> makeMultiAxis(
      const std::vector<std::unique_ptr<IAxis>>& axes);

  std::vector<AxisFactory> m_axisFactories;
};

/// @brief Multi-dimensional binning description of a compile-time dimension
///
/// Adds a statically sized construction and build API on top of
/// @c MultiAxisFactory , mirroring @c IMultiAxis and @c IMultiAxisXD .
///
/// @tparam DIM number of axes (dimension of the grid)
template <std::size_t DIM>
class MultiAxisFactoryXD : public MultiAxisFactory {
 public:
  static_assert(DIM >= 1 && DIM <= 3,
                "MultiAxisFactoryXD supports 1 to 3 dimensions");

  /// Construct from one axis description per dimension
  /// @param axisFactories the axis descriptions, in axis order
  explicit MultiAxisFactoryXD(std::array<AxisFactory, DIM> axisFactories)
      : MultiAxisFactory(std::vector<AxisFactory>(
            std::make_move_iterator(axisFactories.begin()),
            std::make_move_iterator(axisFactories.end()))) {}

  /// Make the type-erased base overloads visible next to the ones below
  using MultiAxisFactory::buildMultiAxis;

  /// Build a multi-axis from fully specified descriptions
  /// @throws std::domain_error if any description is deferred
  /// @return the created multi-axis of dimension @c DIM
  std::unique_ptr<IMultiAxisXD<DIM>> buildMultiAxis() const {
    return downcast(MultiAxisFactory::buildMultiAxis());
  }

  /// Build a multi-axis from deferred descriptions with the consumer supplied
  /// ranges, boundary types and directions
  /// @param options one option set per axis, in axis order
  /// @throws std::domain_error if any description is fully specified
  /// @throws std::invalid_argument if a direction mismatches
  /// @return the created multi-axis of dimension @c DIM
  std::unique_ptr<IMultiAxisXD<DIM>> buildMultiAxis(
      const std::array<AxisFactory::Options, DIM>& options) const {
    return downcast(MultiAxisFactory::buildMultiAxis(options));
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
          "MultiAxisFactoryXD: unexpected multi-axis dimension");
    }
    multiAxis.release();
    return std::unique_ptr<IMultiAxisXD<DIM>>(xd);
  }
};

/// Type alias for a multi-axis description of dimension 1
using MultiAxisFactory1D = MultiAxisFactoryXD<1>;
/// Type alias for a multi-axis description of dimension 2
using MultiAxisFactory2D = MultiAxisFactoryXD<2>;

/// @brief Build the multi-axis of a surface binning against that surface
///
/// A surface spans two local coordinates, so the binning is strictly
/// two-dimensional; use a single bin in one direction to bin along the other
/// only. Both descriptions are matched to the canonical local axis directions
/// of the surface (see @c Surface::localAxes ), either positionally if
/// neither carries a direction, or by direction if both do. Mixing the two is
/// rejected. The axes are returned in canonical direction order.
///
/// @param multiAxisFactory the binning description to build from
/// @param surface the surface to resolve the ranges against
/// @throws std::invalid_argument if the directions cannot be matched or the
///         surface is unsupported
/// @throws std::domain_error if any description is fully specified instead
///         of deferred
/// @return the created multi-axis, in canonical direction order
std::unique_ptr<IMultiAxis2D> resolveMultiAxis(
    const MultiAxisFactory2D& multiAxisFactory, const Surface& surface);

}  // namespace Acts
