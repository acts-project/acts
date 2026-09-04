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
#include <format>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ActsPodioEdm::detail {

/// Encoded index type for local measurement subspace coordinates.
using SubspaceIndex = std::uint8_t;
/// Largest allowed local subspace index value.
constexpr std::uint8_t kMaxSubspaceIndex = 6u;
/// Maximum number of coordinates in an encoded local subspace.
constexpr std::size_t kMaxSubspaceSize = 6u;

/// Encode a list of bound parameter indices into a 32-bit integer.
///
/// This function bit-packs up to `kMaxSubspaceSize` parameter indices into a
/// single 32-bit unsigned integer for storage in EDM4hep format.
///
/// Bit layout:
///   - Bits 0-3:   Number of indices (size)
///   - Bits 4-7:   First index
///   - Bits 8-11:  Second index
///   - Bits 12-15: Third index
///   - (and so on, up to `kMaxSubspaceSize` indices total)
///
/// @param indices Span of parameter indices to encode (max
/// `kMaxSubspaceSize` elements)
/// @return Packed 32-bit unsigned integer containing all indices
inline std::int32_t encodeIndices(std::span<const std::uint8_t> indices) {
  if (indices.size() > kMaxSubspaceSize) {
    throw std::runtime_error(
        std::format("Number of indices exceeds maximum of {} for EDM4hep",
                    kMaxSubspaceSize));
  }

  std::int32_t result = 0;
  std::uint8_t shift = 0;
  result |= (indices.size() << 0);
  shift += 4;

  for (std::uint8_t index : indices) {
    if (index > kMaxSubspaceIndex) {
      throw std::runtime_error(std::format(
          "Index out of range: maximum allowed is {}", kMaxSubspaceIndex));
    }
    result |= (static_cast<std::int32_t>(index) << shift);
    shift += 4;
  }
  return result;
}

/// Decode a 32-bit integer back into a list of bound parameter indices.
///
/// This function unpacks a bit-packed integer (created by encodeIndices) back
/// into the original list of parameter indices. See encodeIndices for the bit
/// layout specification.
///
/// @param type Packed 32-bit unsigned integer containing encoded indices
/// @return Vector of decoded parameter indices (each in
/// `[0, kMaxSubspaceIndex]`)
inline std::vector<SubspaceIndex> decodeIndices(std::int32_t type) {
  std::vector<SubspaceIndex> result;
  std::uint8_t size = type & 0xF;
  if (size > kMaxSubspaceSize) {
    throw std::runtime_error(
        std::format("Number of indices exceeds maximum of {} for EDM4hep",
                    kMaxSubspaceSize));
  }
  result.resize(size);
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = (type >> ((i + 1) * 4)) & 0xF;
    if (result[i] > kMaxSubspaceIndex) {
      throw std::runtime_error(std::format(
          "Index out of range: maximum allowed is {}", kMaxSubspaceIndex));
    }
  }
  return result;
}

/// Find the measurement-vector position of an enum-coded subspace coordinate.
///
/// Searches @p indices for the given enum value (cast to SubspaceIndex) and
/// returns its position. This position is the index into the measurement and
/// covariance vectors stored in a TrackerHitLocal.
///
/// @tparam E  Enum type whose underlying value encodes a subspace index
/// @param indices  Active subspace indices (e.g. from decodeIndices())
/// @param enumVal  The enum value to look up
/// @return  Position of @p enumVal in @p indices
/// @throws std::runtime_error if @p enumVal is not present in @p indices
template <typename E>
  requires std::is_enum_v<E>
inline std::size_t findSubspaceIndex(std::span<const SubspaceIndex> indices,
                                     E enumVal) {
  const auto target = static_cast<SubspaceIndex>(enumVal);
  for (std::size_t i = 0; i < indices.size(); ++i) {
    if (indices[i] == target) {
      return i;
    }
  }
  throw std::runtime_error(
      std::format("Enum value {} not found in subspace indices",
                  static_cast<int>(enumVal)));
}

}  // namespace ActsPodioEdm::detail
