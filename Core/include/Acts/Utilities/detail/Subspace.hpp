// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace Acts::detail {

/// @defgroup subspace Linear subspace definitions
///
/// All types in this group define a linear subspace of a larger vector space
/// with the following properties: the subspace is defined by a set of axis
/// indices in the full space, i.e. there is no need for a non-trivial
/// projection matrix, all subspace axis indices are unique and there are no
/// duplicated ones, and the set of axis indices must be ordered, i.e. the axis
/// order in the subspace is the same as in the full space. The last requirement
/// is not strictly required by the implementation, but is still added to
/// simplify reasoning.
///
/// Only the size of the subspace are defined by the types at compile time.
/// Which axes comprise the subspace is defined at runtime. This allows to use
/// fixed-size computations (selected at compile time) but avoids the
/// combinatorial explosion of also having to handle all possible combination of
/// axes at compile time. This was tried previously and resulted in sizable
/// resource consumption at compile time without runtime benefits.
///
/// All types are intentionally using `std::size_t` as their template values,
/// instead of the more specific index enums, to reduce the number of templates.
/// This is fully compatible as the index enums are required to be convertible
/// to `std::size_t`.
///
/// All types intentionally only define the subspace but not how vectors
/// and matrices are stored to avoid unnecessary coupling between components,
/// i.e here between the pure definition and the storage.
///
/// All types provide `.projectVector(...)` and `.exandVector(...)` methods to
/// convert to/from the subspace. They also provide `.projector()` and
/// `.expander()` methods to create projection and expansion matrices if they
/// are required explicitly. For the specific subspace requirements listed
/// above, the projection and expansion matrices are transpose to each other. In
/// the general case, this does not have to be the case and users are encouraged
/// to use `.projector()` and `.expander()` instead of e.g.
/// `.projector().transpose()`.
///
/// @{

/// Fixed-size subspace representation.
///
/// @tparam kFullSize Size of the full vector space
/// @tparam kSize Size of the subspace
template <std::size_t kFullSize, std::size_t kSize>
class FixedSizeSubspace {
  static_assert(kFullSize <= static_cast<std::size_t>(
                                 std::numeric_limits<std::uint8_t>::max()),
                "Full vector space size is larger than the supported range");
  static_assert(1u <= kSize, "Subspace size must be at least 1");
  static_assert(kSize <= kFullSize,
                "Subspace can only be as large as the full space");

  template <typename source_t>
  using SubspaceVectorFor = Eigen::Matrix<typename source_t::Scalar, kSize, 1>;
  template <typename source_t>
  using FullspaceVectorFor =
      Eigen::Matrix<typename source_t::Scalar, kFullSize, 1>;
  template <typename scalar_t>
  using ProjectionMatrix = Eigen::Matrix<scalar_t, kSize, kFullSize>;
  template <typename scalar_t>
  using ExpansionMatrix = Eigen::Matrix<scalar_t, kFullSize, kSize>;

  // the functionality could also be implemented using a std::bitset where each
  // bit corresponds to an axis in the fullspace and set bits indicate which
  // bits make up the subspace. this would be a more compact representation but
  // complicates the implementation since we can not easily iterate over the
  // indices of the subspace. storing the subspace indices directly requires a
  // bit more memory but is easier to work with. for our typical use cases with
  // n<=8, this still takes only 64bit of memory.
  std::array<std::uint8_t, kSize> m_axes;

 public:
  /// Construct from a container of axis indices.
  ///
  /// @tparam index_t Input index type, must be convertible to std::uint8_t
  /// @param indices Unique, ordered indices
  template <typename index_t>
  constexpr explicit FixedSizeSubspace(
      const std::array<index_t, kSize>& indices) {
    for (std::size_t i = 0u; i < kSize; ++i) {
      assert((indices[i] < kFullSize) &&
             "Axis indices must be within the full space");
      if (0u < i) {
        assert((indices[i - 1u] < indices[i]) &&
               "Axis indices must be unique and ordered");
      }
    }
    for (std::size_t i = 0; i < kSize; ++i) {
      m_axes[i] = static_cast<std::uint8_t>(indices[i]);
    }
  }
  // The subset can not be constructed w/o defining its axis indices.
  FixedSizeSubspace() = delete;
  FixedSizeSubspace(const FixedSizeSubspace&) = default;
  FixedSizeSubspace(FixedSizeSubspace&&) = default;
  FixedSizeSubspace& operator=(const FixedSizeSubspace&) = default;
  FixedSizeSubspace& operator=(FixedSizeSubspace&&) = default;

  /// Size of the subspace.
  static constexpr std::size_t size() { return kSize; }
  /// Size of the full vector space.
  static constexpr std::size_t fullSize() { return kFullSize; }

  /// Access axis index by position.
  ///
  /// @param i Position in the subspace
  /// @return Axis index in the full space
  constexpr std::size_t operator[](std::size_t i) const { return m_axes[i]; }

  /// Axis indices that comprise the subspace.
  ///
  /// The specific container and index type should be considered an
  /// implementation detail. Users should treat the return type as a generic
  /// container whose elements are convertible to `std::size_t`.
  constexpr const std::array<std::uint8_t, kSize>& indices() const {
    return m_axes;
  }

  /// Check if the given axis index in the full space is part of the subspace.
  constexpr bool contains(std::size_t index) const {
    bool isContained = false;
    // always iterate over all elements to avoid branching and hope the compiler
    // can optimise this for us.
    for (std::uint8_t a : m_axes) {
      isContained = isContained || (a == index);
    }
    return isContained;
  }

  /// Project a full vector into the subspace.
  ///
  /// @tparam fullspace_vector_t Vector type in the full space
  /// @param full Vector in the full space
  /// @return Subspace vector w/ just the configured axis components
  ///
  /// @note Always returns a column vector regardless of the input
  template <typename fullspace_vector_t>
  SubspaceVectorFor<fullspace_vector_t> projectVector(
      const Eigen::MatrixBase<fullspace_vector_t>& full) const {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(fullspace_vector_t, kFullSize);

    SubspaceVectorFor<fullspace_vector_t> sub;
    for (std::size_t i = 0u; i < kSize; ++i) {
      sub[i] = full[m_axes[i]];
    }
    return sub;
  }

  /// Expand a subspace vector into the full space.
  ///
  /// @tparam subspace_vector_t Subspace vector type
  /// @param sub Subspace vector
  /// @return Vector in the full space w/ zeros for non-subspace components
  ///
  /// @note Always returns a column vector regardless of the input
  template <typename subspace_vector_t>
  FullspaceVectorFor<subspace_vector_t> expandVector(
      const Eigen::MatrixBase<subspace_vector_t>& sub) const {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(subspace_vector_t, kSize);

    FullspaceVectorFor<subspace_vector_t> full;
    full.setZero();
    for (std::size_t i = 0u; i < kSize; ++i) {
      full[m_axes[i]] = sub[i];
    }
    return full;
  }

  /// Projection matrix that maps from the full space into the subspace.
  ///
  /// @tparam scalar_t Scalar type for the projection matrix
  template <typename scalar_t>
  ProjectionMatrix<scalar_t> projector() const {
    ProjectionMatrix<scalar_t> proj;
    proj.setZero();
    for (std::size_t i = 0u; i < kSize; ++i) {
      proj(i, m_axes[i]) = 1;
    }
    return proj;
  }

  /// Expansion matrix that maps from the subspace into the full space.
  ///
  /// @tparam scalar_t Scalar type of the generated expansion matrix
  template <typename scalar_t>
  ExpansionMatrix<scalar_t> expander() const {
    ExpansionMatrix<scalar_t> expn;
    expn.setZero();
    for (std::size_t i = 0u; i < kSize; ++i) {
      expn(m_axes[i], i) = 1;
    }
    return expn;
  }
};

/// Variable-size subspace representation.
///
/// @tparam kFullSize Size of the full vector space
/// @tparam kSize Size of the subspace
template <std::size_t kFullSize>
class VariableSizeSubspace {
  static_assert(kFullSize <= static_cast<std::size_t>(UINT8_MAX),
                "Full vector space size is larger than the supported range");

  template <typename scalar_t>
  using ProjectionMatrix = Eigen::Matrix<scalar_t, Eigen::Dynamic, kFullSize>;
  template <typename scalar_t>
  using ExpansionMatrix = Eigen::Matrix<scalar_t, kFullSize, Eigen::Dynamic>;

  template <typename scalar_t>
  using FullProjectionMatrix = Eigen::Matrix<scalar_t, kFullSize, kFullSize>;
  template <typename scalar_t>
  using FullExpansionMatrix = Eigen::Matrix<scalar_t, kFullSize, kFullSize>;

  std::size_t m_size{};

  // the functionality could also be implemented using a std::bitset where each
  // bit corresponds to an axis in the fullspace and set bits indicate which
  // bits make up the subspace. this would be a more compact representation but
  // complicates the implementation since we can not easily iterate over the
  // indices of the subspace. storing the subspace indices directly requires a
  // bit more memory but is easier to work with. for our typical use cases with
  // n<=8, this still takes only 64bit of memory.
  std::array<std::uint8_t, kFullSize> m_axes{};

 public:
  /// Construct from a container of axis indices.
  ///
  /// @tparam index_t Input index type, must be convertible to std::uint8_t
  /// @param indices Unique, ordered indices
  template <typename index_t, std::size_t kSize>
  constexpr explicit VariableSizeSubspace(
      const std::array<index_t, kSize>& indices) {
    m_size = kSize;
    for (std::size_t i = 0u; i < kSize; ++i) {
      assert((indices[i] < kFullSize) &&
             "Axis indices must be within the full space");
      if (0u < i) {
        assert((indices[i - 1u] < indices[i]) &&
               "Axis indices must be unique and ordered");
      }
    }
    for (std::size_t i = 0; i < kSize; ++i) {
      m_axes[i] = static_cast<std::uint8_t>(indices[i]);
    }
  }
  // The subset can not be constructed w/o defining its axis indices.
  VariableSizeSubspace() = delete;
  VariableSizeSubspace(const VariableSizeSubspace&) = default;
  VariableSizeSubspace(VariableSizeSubspace&&) = default;
  VariableSizeSubspace& operator=(const VariableSizeSubspace&) = default;
  VariableSizeSubspace& operator=(VariableSizeSubspace&&) = default;

  /// Size of the subspace.
  constexpr std::size_t size() const { return m_size; }
  /// Size of the full vector space.
  static constexpr std::size_t fullSize() { return kFullSize; }

  /// Access axis index by position.
  ///
  /// @param i Position in the subspace
  /// @return Axis index in the full space
  constexpr std::size_t operator[](std::size_t i) const {
    assert(i < m_size);
    return m_axes[i];
  }

  /// Check if the given axis index in the full space is part of the subspace.
  constexpr bool contains(std::size_t index) const {
    bool isContained = false;
    // always iterate over all elements to avoid branching and hope the compiler
    // can optimise this for us.
    for (std::size_t i = 0; i < kFullSize; ++i) {
      isContained = isContained || ((i < m_size) && (m_axes[i] == index));
    }
    return isContained;
  }

  /// Projection matrix that maps from the full space into the subspace.
  ///
  /// @tparam scalar_t Scalar type for the projection matrix
  template <typename scalar_t>
  ProjectionMatrix<scalar_t> projector() const {
    ProjectionMatrix<scalar_t> proj(m_size, kFullSize);
    proj.setZero();
    for (std::size_t i = 0u; i < m_size; ++i) {
      proj(i, m_axes[i]) = 1;
    }
    return proj;
  }

  /// Expansion matrix that maps from the subspace into the full space.
  ///
  /// @tparam scalar_t Scalar type of the generated expansion matrix
  template <typename scalar_t>
  ExpansionMatrix<scalar_t> expander() const {
    ExpansionMatrix<scalar_t> expn(kFullSize, m_size);
    expn.setZero();
    for (std::size_t i = 0u; i < m_size; ++i) {
      expn(m_axes[i], i) = 1;
    }
    return expn;
  }

  /// Projection matrix that maps from the full space into the subspace.
  ///
  /// @tparam scalar_t Scalar type for the projection matrix
  template <typename scalar_t>
  FullProjectionMatrix<scalar_t> fullProjector() const {
    FullProjectionMatrix<scalar_t> proj;
    proj.setZero();
    for (std::size_t i = 0u; i < m_size; ++i) {
      proj(i, m_axes[i]) = 1;
    }
    return proj;
  }

  /// Expansion matrix that maps from the subspace into the full space.
  ///
  /// @tparam scalar_t Scalar type of the generated expansion matrix
  template <typename scalar_t>
  FullExpansionMatrix<scalar_t> fullExpander() const {
    FullExpansionMatrix<scalar_t> expn;
    expn.setZero();
    for (std::size_t i = 0u; i < m_size; ++i) {
      expn(m_axes[i], i) = 1;
    }
    return expn;
  }

  std::uint64_t projectorBits() const {
    std::uint64_t result = 0;

    for (std::size_t i = 0; i < m_size; ++i) {
      for (std::size_t j = 0; j < kFullSize; ++j) {
        // the bit order is defined in `Acts/Utilities/AlgebraHelpers.hpp`
        // in `matrixToBitset`
        std::size_t index = m_size * kFullSize - 1 - (i + j * m_size);
        if (m_axes[i] == j) {
          result |= (1ull << index);
        }
      }
    }

    return result;
  }

  std::uint64_t fullProjectorBits() const {
    std::uint64_t result = 0;

    for (std::size_t i = 0; i < kFullSize; ++i) {
      for (std::size_t j = 0; j < kFullSize; ++j) {
        // the bit order is defined in `Acts/Utilities/AlgebraHelpers.hpp`
        // in `matrixToBitset`
        std::size_t index = kFullSize * kFullSize - 1 - (i + j * kFullSize);
        if (i < m_size && m_axes[i] == j) {
          result |= (1ull << index);
        }
      }
    }

    return result;
  }
};

/// @}

}  // namespace Acts::detail
