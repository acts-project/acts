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
#include <iosfwd>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace Acts {

/// @brief Description of a multi-dimensional binning that produces the axes
/// of a grid
///
/// This class bundles one @c AxisFactory per grid dimension. Analogous to the
/// @c IMultiAxis / @c IMultiAxisXD split, the dimension is only known at
/// runtime through this base class, while the dimension-aware API is exposed
/// by the derived @c MultiAxisFactoryXD template. Unlike @c IMultiAxis the
/// class is value-semantic: the axis descriptions are stored in this base
/// class, so it can be held by value e.g. in configuration objects.
///
/// The resolution API mirrors @c AxisFactory : fully specified descriptions
/// resolve without input, deferred descriptions require one @c AxisResolution
/// per axis, supplied positionally in storage order.
class MultiAxisFactory {
 public:
  /// Construct from one axis description per dimension
  /// @param axisFactories the axis descriptions, in axis order
  /// @throws std::invalid_argument if no axis description is given
  explicit MultiAxisFactory(std::vector<AxisFactory> axisFactories);

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
  /// consumer supplied resolutions to produce axes
  /// @return true if any axis description is deferred
  bool isDeferred() const;

  /// Produce the axes from fully specified descriptions
  /// @throws std::domain_error if any description is deferred
  /// @return the created axes, in axis order
  std::vector<std::unique_ptr<IAxis>> toAxes() const;

  /// Produce the axes from deferred descriptions with the consumer supplied
  /// ranges and boundary types
  /// @param resolutions one range and boundary type per axis, in axis order
  /// @param directions optionally one caller expected direction per axis, in
  ///        axis order; if given, each has to agree with the corresponding
  ///        stored direction
  /// @throws std::domain_error if any description is fully specified
  /// @throws std::invalid_argument if the number of resolutions or directions
  ///         does not match the number of axes, or a direction mismatches
  /// @return the created axes, in axis order
  std::vector<std::unique_ptr<IAxis>> toAxes(
      std::span<const AxisResolution> resolutions,
      std::span<const AxisDirection> directions = {}) const;

  /// Produce a multi-axis from fully specified descriptions
  /// @throws std::domain_error if any description is deferred or the
  ///         dimension exceeds the supported maximum of 3
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> toMultiAxis() const;

  /// Produce a multi-axis from deferred descriptions with the consumer
  /// supplied ranges and boundary types
  /// @param resolutions one range and boundary type per axis, in axis order
  /// @param directions optionally one caller expected direction per axis, in
  ///        axis order; if given, each has to agree with the corresponding
  ///        stored direction
  /// @throws std::domain_error if any description is fully specified or the
  ///         dimension exceeds the supported maximum of 3
  /// @throws std::invalid_argument if the number of resolutions or directions
  ///         does not match the number of axes, or a direction mismatches
  /// @return the created multi-axis
  std::unique_ptr<IMultiAxis> toMultiAxis(
      std::span<const AxisResolution> resolutions,
      std::span<const AxisDirection> directions = {}) const;

  /// Get a string representation of this description
  /// @return the string representation
  std::string toString() const;

  /// Check if two descriptions are equal
  /// @param lhs first description
  /// @param rhs second description
  /// @return true if all axis descriptions are equal
  friend bool operator==(const MultiAxisFactory& lhs,
                         const MultiAxisFactory& rhs) {
    return lhs.m_axisFactories == rhs.m_axisFactories;
  }

  /// Output stream operator
  /// @param os output stream
  /// @param multiAxisFactory the description to be printed
  /// @return the output stream
  friend std::ostream& operator<<(std::ostream& os,
                                  const MultiAxisFactory& multiAxisFactory);

 private:
  /// Create a multi-axis from already produced axes
  /// @param axes the axes, in axis order
  /// @return the created multi-axis
  static std::unique_ptr<IMultiAxis> makeMultiAxis(
      const std::vector<std::unique_ptr<IAxis>>& axes);

  std::vector<AxisFactory> m_axisFactories;
};

/// @brief Multi-dimensional binning description of a fixed, compile-time
/// dimension
///
/// On top of the runtime-dimension @c MultiAxisFactory API this adds a
/// statically sized construction and resolution API, mirroring the relation
/// between @c IMultiAxis and @c IMultiAxisXD.
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

  /// Produce a multi-axis from fully specified descriptions
  /// @throws std::domain_error if any description is deferred
  /// @return the created multi-axis of dimension @c DIM
  std::unique_ptr<IMultiAxisXD<DIM>> toMultiAxis() const {
    return downcast(MultiAxisFactory::toMultiAxis());
  }

  /// Produce a multi-axis from deferred descriptions with the consumer
  /// supplied ranges and boundary types
  /// @param resolutions one range and boundary type per axis, in axis order
  /// @param directions optionally one caller expected direction per axis, in
  ///        axis order; if given, each has to agree with the corresponding
  ///        stored direction
  /// @throws std::domain_error if any description is fully specified
  /// @throws std::invalid_argument if a direction mismatches
  /// @return the created multi-axis of dimension @c DIM
  std::unique_ptr<IMultiAxisXD<DIM>> toMultiAxis(
      const std::array<AxisResolution, DIM>& resolutions,
      std::span<const AxisDirection> directions = {}) const {
    return downcast(MultiAxisFactory::toMultiAxis(resolutions, directions));
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

/// Output stream operator
/// @param os output stream
/// @param multiAxisFactory the description to be printed
/// @return the output stream
std::ostream& operator<<(std::ostream& os,
                         const MultiAxisFactory& multiAxisFactory);

/// Type alias for a multi-axis description of dimension 1
using MultiAxisFactory1D = MultiAxisFactoryXD<1>;
/// Type alias for a multi-axis description of dimension 2
using MultiAxisFactory2D = MultiAxisFactoryXD<2>;

}  // namespace Acts
