// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <cstddef>
#include <ranges>

#include <boost/container/static_vector.hpp>

namespace Acts {

/// @brief Check subspace indices for consistency
///
/// Indices must be unique and within the full size of the subspace
///
/// @tparam index_range_t the type of the container of indices
///
/// @param indexRange the range of indices
/// @param fullSize the full size of the subspace
/// @param subspaceSize the size of the subspace
///
/// @return true if the indices are consistent
template <std::ranges::sized_range index_range_t>
inline static bool checkSubspaceIndices(const index_range_t& indexRange,
                                        std::size_t fullSize,
                                        std::size_t subspaceSize) {
  if (subspaceSize > fullSize) {
    return false;
  }
  if (static_cast<std::size_t>(indexRange.size()) != subspaceSize) {
    return false;
  }
  for (auto it = indexRange.begin(); it != indexRange.end();) {
    auto index = *it;
    if (index >= fullSize) {
      return false;
    }
    ++it;
    if (std::find(it, indexRange.end(), index) != indexRange.end()) {
      return false;
    }
  }
  return true;
}

/// Serialize subspace indices to a single 64 bit integer
///
/// @tparam FullSize the full size of the subspace
///
/// @param indices the subspace indices
///
/// @return the serialized subspace indices
template <std::size_t FullSize>
inline static SerializedSubspaceIndices serializeSubspaceIndices(
    const SubspaceIndices<FullSize>& indices)
  requires(FullSize <= 8)
{
  SerializedSubspaceIndices result = 0;
  for (std::size_t i = 0; i < FullSize; ++i) {
    result |= static_cast<SerializedSubspaceIndices>(indices[i] & 0xFF)
              << (i * 8);
  }
  return result;
}

/// Deserialize subspace indices from a single 64 bit integer
///
/// @tparam FullSize the full size of the subspace
///
/// @param serialized the serialized subspace indices
///
/// @return the subspace indices
template <std::size_t FullSize>
inline static SubspaceIndices<FullSize> deserializeSubspaceIndices(
    SerializedSubspaceIndices serialized)
  requires(FullSize <= 8)
{
  SubspaceIndices<FullSize> result;
  for (std::size_t i = 0; i < FullSize; ++i) {
    result[i] = static_cast<std::uint8_t>((serialized >> (i * 8)) & 0xFF);
  }
  return result;
}

namespace detail {

/// Helper base class for subspace operations which provides common
/// functionality for fixed and variable subspace helpers
///
/// @tparam Derived the derived type
/// @tparam FullSize the full size of the subspace
template <typename Derived, std::size_t FullSize>
class SubspaceHelperBase {
 public:
  static constexpr std::size_t kFullSize = FullSize;

  using FullSquareMatrix = SquareMatrix<kFullSize>;

  using size_type = std::size_t;

  bool empty() const { return self().empty(); }
  std::size_t size() const { return self().size(); }

  auto operator[](std::size_t i) const { return self()[i]; }

  auto begin() const { return self().begin(); }
  auto end() const { return self().end(); }

  bool contains(std::uint8_t index) const {
    return std::find(begin(), end(), index) != end();
  }
  std::size_t indexOf(std::uint8_t index) const {
    auto it = std::find(begin(), end(), index);
    return it != end() ? std::distance(begin(), it) : kFullSize;
  }

  template <typename EigenDerived>
  Vector<kFullSize> expandVector(
      const Eigen::DenseBase<EigenDerived>& vector) const {
    Vector<kFullSize> result = Vector<kFullSize>::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(index) = vector(i);
    }
    return result;
  }

  template <typename EigenDerived>
  FullSquareMatrix expandMatrix(
      const Eigen::DenseBase<EigenDerived>& matrix) const {
    FullSquareMatrix result = FullSquareMatrix::Zero();
    for (auto [i, indexI] : enumerate(*this)) {
      for (auto [j, indexJ] : enumerate(*this)) {
        result(indexI, indexJ) = matrix(i, j);
      }
    }
    return result;
  }

  FullSquareMatrix fullProjector() const {
    FullSquareMatrix result = FullSquareMatrix::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(i, index) = 1;
    }
    return result;
  }

  FullSquareMatrix fullExpander() const {
    FullSquareMatrix result = FullSquareMatrix::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(index, i) = 1;
    }
    return result;
  }

  ProjectorBitset projectorBitset() const
    requires(kFullSize <= 8)
  {
    return matrixToBitset(fullProjector()).to_ullong();
  }

 private:
  const Derived& self() const { return static_cast<const Derived&>(*this); }
};

}  // namespace detail

/// Helper class for variable subspace operations
///
/// @tparam FullSize the full size of the subspace
/// @tparam index_t the index type
template <std::size_t FullSize, typename index_t = std::uint8_t>
class VariableSubspaceHelper
    : public detail::SubspaceHelperBase<
          VariableSubspaceHelper<FullSize, index_t>, FullSize> {
 public:
  /// Size of the full parameter space
  static constexpr std::size_t kFullSize = FullSize;

  /// Type alias for index type
  using IndexType = index_t;
  /// Type alias for container holding subspace indices
  using Container = boost::container::static_vector<IndexType, FullSize>;

  /// Construct a variable subspace helper with specified indices
  ///
  /// @tparam other_index_range_t Type of the index range
  /// @param indices Range of indices defining the subspace
  template <std::ranges::sized_range other_index_range_t>
  explicit VariableSubspaceHelper(const other_index_range_t& indices) {
    assert(checkSubspaceIndices(indices, kFullSize, indices.size()) &&
           "Invalid indices");
    m_indices.resize(indices.size());
    std::transform(indices.begin(), indices.end(), m_indices.begin(),
                   [](auto index) { return static_cast<IndexType>(index); });
  }

  /// Check if the subspace is empty
  /// @return True if the subspace contains no indices
  bool empty() const { return m_indices.empty(); }
  /// Get the size of the subspace
  /// @return Number of indices in the subspace
  std::size_t size() const { return m_indices.size(); }
  /// Get the container of subspace indices
  /// @return Reference to the container holding the indices
  const Container& indices() const { return m_indices; }

  /// Access subspace index at position i
  /// @param i Position index
  /// @return The subspace index at position i
  IndexType operator[](std::size_t i) const { return m_indices[i]; }

  /// Get iterator to beginning of indices
  /// @return Iterator to the first index
  auto begin() const { return m_indices.begin(); }
  /// Get iterator to end of indices
  /// @return Iterator past the last index
  auto end() const { return m_indices.end(); }

 private:
  Container m_indices;
};

/// Helper class for fixed subspace operations
///
/// @tparam FullSize the full size of the subspace
/// @tparam SubspaceSize the size of the subspace
/// @tparam index_t the index type
template <std::size_t FullSize, std::size_t SubspaceSize,
          typename index_t = std::uint8_t>
class FixedSubspaceHelper
    : public detail::SubspaceHelperBase<
          FixedSubspaceHelper<FullSize, SubspaceSize, index_t>, FullSize> {
 public:
  /// Size of the full parameter space
  static constexpr std::size_t kFullSize = FullSize;
  /// Size of the subspace parameter space
  static constexpr std::size_t kSubspaceSize = SubspaceSize;
  static_assert(kSubspaceSize <= kFullSize, "Invalid subspace size");

  /// Type alias for projection matrix (subspace x full)
  using Projector = Matrix<kSubspaceSize, kFullSize>;
  /// Type alias for expansion matrix (full x subspace)
  using Expander = Matrix<kFullSize, kSubspaceSize>;
  /// Type alias for subspace vector
  using VectorD = Vector<kSubspaceSize>;
  /// Type alias for subspace square matrix
  using SquareMatrixD = SquareMatrix<kSubspaceSize>;
  /// Type alias for left application result matrix
  template <std::size_t K>
  using ApplyLeftResult = Matrix<kSubspaceSize, kSubspaceSize>;
  /// Type alias for right application result matrix
  template <std::size_t N>
  using ApplyRightResult = Matrix<kSubspaceSize, kSubspaceSize>;

  /// Type alias for index type used to specify subspace indices
  using IndexType = index_t;
  /// Type alias for container storing subspace indices
  using Container = std::array<IndexType, kSubspaceSize>;

  /// Type alias for size type
  using size_type = std::size_t;

  /// Construct a fixed subspace helper with specified indices
  ///
  /// @tparam other_index_range_t Type of the index range
  /// @param indices Range of indices defining the subspace, must match SubspaceSize
  template <std::ranges::sized_range other_index_range_t>
  explicit FixedSubspaceHelper(const other_index_range_t& indices) {
    assert(checkSubspaceIndices(indices, kFullSize, kSubspaceSize) &&
           "Invalid indices");
    std::transform(indices.begin(), indices.end(), m_indices.begin(),
                   [](auto index) { return static_cast<IndexType>(index); });
  }

  /// Check if the subspace is empty
  /// @return True if the subspace contains no indices (always false for fixed subspaces)
  bool empty() const { return m_indices.empty(); }
  /// Get the size of the subspace
  /// @return Number of indices in the subspace (always SubspaceSize)
  std::size_t size() const { return m_indices.size(); }
  /// Get the container of subspace indices
  /// @return Reference to the array holding the indices
  const Container& indices() const { return m_indices; }

  /// Access subspace index at position i
  /// @param i Position index
  /// @return The subspace index at position i
  IndexType operator[](std::uint32_t i) const { return m_indices[i]; }

  /// Get iterator to beginning of indices
  /// @return Iterator to the first index
  auto begin() const { return m_indices.begin(); }
  /// Get iterator to end of indices
  /// @return Iterator past the last index
  auto end() const { return m_indices.end(); }

  /// Create a projection matrix from full space to subspace
  /// @return Matrix that projects from full space to subspace dimensions
  Projector projector() const {
    Projector result = Projector::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(i, index) = 1;
    }
    return result;
  }

  /// Create an expansion matrix from subspace to full space
  /// @return Matrix that expands from subspace to full space dimensions
  Expander expander() const {
    Expander result = Expander::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(index, i) = 1;
    }
    return result;
  }

  /// Apply subspace projection on the left side of a matrix
  /// @tparam K Number of columns in the input matrix
  /// @tparam Derived Eigen matrix type
  /// @param matrix Input matrix with rows matching full space size
  /// @return Matrix with subspace rows and K columns
  template <std::size_t K, typename Derived>
  ApplyLeftResult<K> applyLeftOf(
      const Eigen::DenseBase<Derived>& matrix) const {
    assert(matrix.rows() == kFullSize && "Invalid matrix size");
    ApplyLeftResult<K> result = ApplyLeftResult<K>::Zero();
    for (auto [i, indexI] : enumerate(*this)) {
      for (std::size_t j = 0; j < K; ++j) {
        result(i, j) = matrix(indexI, j);
      }
    }
    return result;
  }

  /// Apply subspace projection on the right side of a matrix
  /// @tparam N Number of rows in the input matrix
  /// @tparam Derived Eigen matrix type
  /// @param matrix Input matrix with columns matching subspace size
  /// @return Matrix with N rows and subspace columns
  template <std::size_t N, typename Derived>
  ApplyRightResult<N> applyRightOf(
      const Eigen::DenseBase<Derived>& matrix) const {
    assert(matrix.cols() == kSubspaceSize && "Invalid matrix size");
    ApplyRightResult<N> result = ApplyRightResult<N>::Zero();
    for (std::size_t i = 0; i < N; ++i) {
      for (auto [j, indexJ] : enumerate(*this)) {
        result(i, j) = matrix(i, indexJ);
      }
    }
    return result;
  }

  /// Project a full-space vector to subspace
  /// @tparam Derived Eigen vector type
  /// @param fullVector Input vector with full space dimensions
  /// @return Vector projected to subspace dimensions
  template <typename Derived>
  VectorD projectVector(const Eigen::DenseBase<Derived>& fullVector) const {
    assert(fullVector.size() == kFullSize && "Invalid full vector size");
    VectorD result = VectorD::Zero();
    for (auto [i, index] : enumerate(*this)) {
      result(i) = fullVector(index);
    }
    return result;
  }

  /// Project a full-space square matrix to subspace
  /// @tparam Derived Eigen matrix type
  /// @param fullMatrix Input square matrix with full space dimensions
  /// @return Square matrix projected to subspace dimensions
  template <typename Derived>
  SquareMatrixD projectMatrix(
      const Eigen::DenseBase<Derived>& fullMatrix) const {
    assert(fullMatrix.rows() == kFullSize && fullMatrix.cols() == kFullSize &&
           "Invalid full matrix size");
    SquareMatrixD result = SquareMatrixD::Zero();
    for (auto [i, indexI] : enumerate(*this)) {
      for (auto [j, indexJ] : enumerate(*this)) {
        result(i, j) = fullMatrix(indexI, indexJ);
      }
    }
    return result;
  }

 private:
  Container m_indices{};
};

/// @brief Helper type for fixed-size bound parameter subspaces
/// @details Provides utilities for working with fixed-size subspaces of bound track parameters
template <std::size_t SubspaceSize>
using FixedBoundSubspaceHelper =
    FixedSubspaceHelper<Acts::eBoundSize, SubspaceSize, std::uint8_t>;

/// @typedef VariableBoundSubspaceHelper
/// Helper type for variable-size bound parameter subspaces.
/// Provides utilities for working with variable-size subspaces of bound track
/// parameters.
using VariableBoundSubspaceHelper =
    VariableSubspaceHelper<Acts::eBoundSize, std::uint8_t>;

/// Convert a projector to subspace indices
///
/// @tparam kFullSize the full size of the subspace
/// @tparam Derived the derived type
///
/// @param projector the projector
///
/// @return the subspace indices
template <std::size_t kFullSize, typename Derived>
SubspaceIndices<kFullSize> projectorToSubspaceIndices(
    const Eigen::DenseBase<Derived>& projector) {
  auto rows = static_cast<std::size_t>(projector.rows());
  auto cols = static_cast<std::size_t>(projector.cols());
  assert(cols == kFullSize && rows <= kFullSize && "Invalid projector size");
  SubspaceIndices<kFullSize> result;
  result.fill(kFullSize);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      assert((projector(i, j) == 0 || projector(i, j) == 1) &&
             "Invalid projector value");
      if (projector(i, j) == 1) {
        result[i] = j;
      }
    }
  }
  return result;
}

}  // namespace Acts
