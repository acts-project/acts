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

namespace Acts {
namespace detail {

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
/// All types are intentionally using `size_t` as their template values, instead
/// of the more specific index enums, to reduce the number of templates. This is
/// fully compatible as the index enums are required to be convertible to
/// `size_t`.
///
/// All types intentionally only define the subspace but not how vectors
/// and matrices are stored to avoid unnecessary coupling between components,
/// i.e here between the pure definition and the storage.
///
/// All types provide `.projectVector(...)` and `.exandVector(...)` methods to
/// convert to/from the subspace. They also provide `.projector()` and
/// `.expander()` methods to create projection and expansion matrices if they
/// are required explicitely. For the specific subspace requirements listed
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
template <size_t kFullSize, size_t kSize>
class FixedSizeSubspace {
  static_assert(kFullSize <= static_cast<size_t>(UINT8_MAX),
                "Full vector space size is larger than the supported range");
  static_assert(1u <= kSize, "Subspace size must be at least 1");
  static_assert(kSize <= kFullSize,
                "Subspace can only be as large as the full space");

  template <typename source_t>
  using SubspaceVectorFor = ActsVector<typename source_t::Scalar, kSize>;
  template <typename source_t>
  using FullspaceVectorFor = ActsVector<typename source_t::Scalar, kFullSize>;
  template <typename scalar_t>
  using ProjectionMatrix = ActsMatrix<scalar_t, kSize, kFullSize>;
  template <typename scalar_t>
  using ExpansionMatrix = ActsMatrix<scalar_t, kFullSize, kSize>;

  // the functionality could also be implemented using a std::bitset where each
  // bit corresponds to an axis in the fullspace and set bits indicate which
  // bits make up the subspace. this would be a more compact representation but
  // complicates the implementation since we can not easily iterate over the
  // indices of the subspace. storing the subspace indices directly requires a
  // bit more memory but is easier to work with. for our typical use cases with
  // n<=8, this still takes only 64bit of memory.
  std::array<uint8_t, kSize> m_axes;

 public:
  /// Construct from a container of axis indices.
  ///
  /// @tparam index_container_t Forward-iterable Container for indices
  /// @param indices Ordered, unique axis indices, must contain kSize entries
  template <typename index_container_t>
  constexpr FixedSizeSubspace(const index_container_t& indices) {
    // this way of iteration imposes least constraints on the input container
    auto idx = std::begin(indices);
    auto end = std::end(indices);
    size_t istore = 0u;
    // limit iteration both by input and storage size
    for (; (idx != end) and (istore < kSize); ++idx, ++istore) {
      // verify inputs
      assert((*idx < kFullSize) and
             "Axis indices must be within the full space");
      if (0u < istore) {
        assert((m_axes[istore - 1u] < *idx) and
               "Axis indices must be strongly ordered and unique");
      }
      // assign stored axis
      m_axes[istore] = static_cast<uint8_t>(*idx);
    }
    assert((istore == kSize) and "Not enough input indices");
    assert((idx == end) and "Too many input indices");
    // if there are missing input indices, fill the axes with invalid values
    // this is a defensive measure for non-debug builds to provoke segfaults
    for (; istore < kSize; ++istore) {
      m_axes[istore] = UINT8_MAX;
    }
  }
  // The subset can not be constructed w/o defining its axis indices.
  FixedSizeSubspace() = delete;
  FixedSizeSubspace(const FixedSizeSubspace&) = default;
  FixedSizeSubspace(FixedSizeSubspace&&) = default;
  FixedSizeSubspace& operator=(const FixedSizeSubspace&) = default;
  FixedSizeSubspace& operator=(FixedSizeSubspace&&) = default;

  /// Size of the subspace.
  static constexpr size_t size() { return kSize; }
  /// Size of the full vector space.
  static constexpr size_t fullSize() { return kFullSize; }

  /// Axis indices that comprise the subspace.
  ///
  /// The specific container and index type should be considered an
  /// implementation detail. Users should treat the return type as a generic
  /// container whose elements are convertible to `size_t`.
  constexpr const std::array<uint8_t, kSize>& indices() const { return m_axes; }

  /// Check if the given axis index in the full space is part of the subspace.
  constexpr bool contains(size_t index) const {
    bool isContained = false;
    // always iterate over all elements to avoid branching and hope the compiler
    // can optimized this for us.
    for (auto a : m_axes) {
      isContained = (isContained or (a == index));
    }
    return isContained;
  }

  /// Project a full vector into the subspace.
  ///
  /// @tparam fullspace_vector_t Veector type in the full space
  /// @param full Vector in the full space
  /// @return Subspace vector w/ just the configured axis components
  ///
  /// @note Always returns a column vector regardless of the input
  template <typename fullspace_vector_t>
  auto projectVector(const Eigen::MatrixBase<fullspace_vector_t>& full) const
      -> SubspaceVectorFor<fullspace_vector_t> {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(fullspace_vector_t, kFullSize);

    SubspaceVectorFor<fullspace_vector_t> sub;
    for (auto i = 0u; i < kSize; ++i) {
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
  auto expandVector(const Eigen::MatrixBase<subspace_vector_t>& sub) const
      -> FullspaceVectorFor<subspace_vector_t> {
    EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(subspace_vector_t, kSize);

    FullspaceVectorFor<subspace_vector_t> full;
    full.setZero();
    for (auto i = 0u; i < kSize; ++i) {
      full[m_axes[i]] = sub[i];
    }
    return full;
  }

  /// Projection matrix that maps from the full space into the subspace.
  ///
  /// @tparam scalar_t Scalar type for the projection matrix
  template <typename scalar_t>
  auto projector() const -> ProjectionMatrix<scalar_t> {
    ProjectionMatrix<scalar_t> proj;
    proj.setZero();
    for (auto i = 0u; i < kSize; ++i) {
      proj(i, m_axes[i]) = 1;
    }
    return proj;
  }

  /// Expansion matrix that maps from the subspace into the full space.
  ///
  /// @tparam scalar_t Scalar type of the generated expansion matrix
  template <typename scalar_t>
  auto expander() const -> ExpansionMatrix<scalar_t> {
    ExpansionMatrix<scalar_t> expn;
    expn.setZero();
    for (auto i = 0u; i < kSize; ++i) {
      expn(m_axes[i], i) = 1;
    }
    return expn;
  }
};

/// @}

}  // namespace detail
}  // namespace Acts
