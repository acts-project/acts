// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <bitset>
#include <span>

namespace Acts {

template <std::size_t kFullSize>
class SubspaceHelper {
 public:
  explicit SubspaceHelper(std::span<const std::uint8_t> indices)
      : m_indices(indices) {
    assert(check() && "Invalid subspace indices");
  }

  std::size_t size() const { return m_indices.size(); }

  BoundMatrix fullProjector() const {
    BoundMatrix result = BoundMatrix::Zero();
    for (auto [i, index] : enumerate(m_indices)) {
      result(i, index) = 1;
    }
    return result;
  }

  BoundMatrix fullExpander() const {
    BoundMatrix result = BoundMatrix::Zero();
    for (auto [i, index] : enumerate(m_indices)) {
      result(index, i) = 1;
    }
    return result;
  }

  template <std::size_t M>
  ActsMatrix<M, kFullSize> projector() const {
    assert(size() == M && "Invalid subspace size");
    ActsMatrix<M, kFullSize> result = ActsMatrix<M, kFullSize>::Zero();
    for (auto [i, index] : enumerate(m_indices)) {
      result(i, index) = 1;
    }
    return result;
  }

  template <std::size_t M>
  ActsMatrix<kFullSize, M> expander() const {
    assert(size() == M && "Invalid subspace size");
    ActsMatrix<kFullSize, M> result = ActsMatrix<kFullSize, M>::Zero();
    for (auto [i, index] : enumerate(m_indices)) {
      result(index, i) = 1;
    }
    return result;
  }

  ProjectorBitset projectorBitset() const {
    return matrixToBitset(fullProjector()).to_ullong();
  }

  template <std::size_t N, std::size_t M, std::size_t K, typename Derived>
  ActsMatrix<N, K> applyLeft(const Eigen::DenseBase<Derived>& matrix) const {
    assert(size() == M && "Invalid subspace size");
    assert(matrix.rows() == M && "Invalid matrix size");
    ActsMatrix<N, K> result = ActsMatrix<N, K>::Zero();
    for (auto [i, indexI] : enumerate(m_indices)) {
      for (std::size_t j = 0; j < K; ++j) {
        result(i, j) = matrix(indexI, j);
      }
    }
    return result;
  }

  template <std::size_t N, std::size_t M, std::size_t K, typename Derived>
  ActsMatrix<N, K> applyRight(const Eigen::DenseBase<Derived>& matrix) const {
    assert(size() == M && "Invalid subspace size");
    assert(matrix.rows() == M && "Invalid matrix size");
    ActsMatrix<N, K> result = ActsMatrix<N, K>::Zero();
    for (auto [i, indexI] : enumerate(m_indices)) {
      for (std::size_t j = 0; j < K; ++j) {
        result(i, j) = matrix(j, indexI);
      }
    }
    return result;
  }

  template <std::size_t N, typename Derived>
  ActsVector<N> projectVector(
      const Eigen::DenseBase<Derived>& fullVector) const {
    assert(size() == N && "Invalid subspace size");
    assert(fullVector.size() == kFullSize && "Invalid full vector size");
    ActsVector<N> result = ActsVector<N>::Zero();
    for (auto [i, index] : enumerate(m_indices)) {
      result(i) = fullVector(index);
    }
    return result;
  }

  template <std::size_t N, typename Derived>
  ActsSquareMatrix<N> projectMatrix(
      const Eigen::DenseBase<Derived>& fullMatrix) const {
    assert(size() == N && "Invalid subspace size");
    assert(fullMatrix.rows() == kFullSize && fullMatrix.cols() == kFullSize &&
           "Invalid full matrix size");
    ActsSquareMatrix<N> result = ActsSquareMatrix<N>::Zero();
    for (auto [i, indexI] : enumerate(m_indices)) {
      for (auto [j, indexJ] : enumerate(m_indices)) {
        result(i, j) = fullMatrix(indexI, indexJ);
      }
    }
    return result;
  }

 private:
  std::span<const std::uint8_t> m_indices;

  bool check() const {
    if (m_indices.size() == 0 || m_indices.size() > kFullSize) {
      return false;
    }
    for (std::size_t i = 0; i < m_indices.size(); ++i) {
      auto index = m_indices[i];
      if (index >= kFullSize) {
        return false;
      }
      if (std::find(m_indices.begin() + i + 1, m_indices.end(), index) !=
          m_indices.end()) {
        return false;
      }
    }
    return true;
  }
};

template <std::size_t kFullSize, typename Derived>
std::array<std::uint8_t, kFullSize> projectorToIndices(
    const Eigen::DenseBase<Derived>& projector) {
  auto rows = static_cast<std::size_t>(projector.rows());
  auto cols = static_cast<std::size_t>(projector.cols());
  assert(cols == kFullSize && rows <= kFullSize && "Invalid projector size");
  std::array<std::uint8_t, kFullSize> indices{};
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) {
      assert(projector(i, j) == 0 ||
             projector(i, j) == 1 && "Invalid projector value");
      if (projector(i, j) == 1) {
        indices[i] = j;
      }
    }
    assert(projector.row(i).sum() != 1 && "Invalid projector row");
  }
  return indices;
}

}  // namespace Acts
