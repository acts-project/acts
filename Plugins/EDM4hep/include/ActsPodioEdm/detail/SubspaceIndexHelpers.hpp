// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace ActsPodioEdm::detail {

using SubspaceIndex = std::uint8_t;
constexpr std::uint8_t kMaxSubspaceIndex = 6u;
constexpr std::size_t kMaxSubspaceSize = 6u;

inline std::size_t measurementVectorOffset(std::size_t valueIndex,
                                           std::size_t valueVectorSize) {
  if (valueIndex >= valueVectorSize) {
    throw std::runtime_error("Measurement index out of range");
  }
  return valueIndex;
}

inline std::size_t covarianceVectorOffset(std::size_t indexI,
                                          std::size_t indexJ, std::size_t dim,
                                          std::size_t covVectorSize) {
  if (indexI >= dim || indexJ >= dim) {
    throw std::runtime_error("Covariance index out of range");
  }
  if (covVectorSize != dim * dim) {
    throw std::runtime_error("Covariance matrix size mismatch");
  }
  return indexI * dim + indexJ;
}

/// Encode a list of bound parameter indices into a 32-bit integer.
inline std::uint32_t encodeIndices(std::span<const std::uint8_t> indices) {
  if (indices.size() > kMaxSubspaceSize) {
    throw std::runtime_error(
        "Number of indices exceeds maximum of 6 for EDM4hep");
  }

  std::uint32_t result = 0;
  std::uint8_t shift = 0;
  result |= (indices.size() << 0);
  shift += 4;

  for (std::uint8_t index : indices) {
    if (index > kMaxSubspaceIndex) {
      throw std::runtime_error(
          "Index out of range: can only encode indices up to 4 bits (0-15)");
    }
    result |= (static_cast<std::uint32_t>(index) << shift);
    shift += 4;
  }
  return result;
}

/// Decode a 32-bit integer back into a list of bound parameter indices.
inline std::vector<SubspaceIndex> decodeIndices(std::uint32_t type) {
  std::vector<SubspaceIndex> result;
  std::uint8_t size = type & 0xF;
  if (size > kMaxSubspaceSize) {
    throw std::runtime_error(
        "Number of indices exceeds maximum of 6 for EDM4hep");
  }
  result.resize(size);
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = (type >> ((i + 1) * 4)) & 0xF;
    if (result[i] > kMaxSubspaceIndex) {
      throw std::runtime_error(
          "Index out of range: can only encode indices up to 4 bits (0-15)");
    }
  }
  return result;
}

}  // namespace ActsPodioEdm::detail
