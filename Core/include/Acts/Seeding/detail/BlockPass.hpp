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

namespace Acts::detail {

/// Evaluates a cut for each of @p n candidates and stores whether each
/// passed. The candidates are processed in chunks of @p kChunk, then the
/// remainder.
///
/// @tparam kChunk the number of candidates per chunk
/// @param n the number of candidates
/// @param passes the cut, called with the candidate's position, returning 1
///   if it passed and 0 if not
/// @param flags receives each candidate's result, 1 or 0
/// @return the number of candidates that passed
template <std::size_t kChunk, typename Passes>
std::uint32_t evaluateCut(std::size_t n, const Passes& passes,
                          std::uint32_t* flags) {
  std::uint32_t nPassed = 0;
  std::size_t c = 0;
  for (; c + kChunk <= n; c += kChunk) {
    for (std::size_t l = 0; l < kChunk; ++l) {
      flags[c + l] = passes(c + l);
      nPassed += flags[c + l];
    }
  }
  for (std::size_t i = c; i < n; ++i) {
    flags[i] = passes(i);
    nPassed += flags[i];
  }
  return nPassed;
}

/// Collects, in order, the values of the candidates whose flag is set. Each
/// candidate's value is written to the next free position, which advances
/// only if its flag is set.
///
/// @param n the number of candidates
/// @param flags whether each candidate is kept, 1 or 0
/// @param valueOf the value collected, called with the candidate's position
/// @param out receives the values of the kept candidates; room for @p n
/// @return the number of candidates kept
template <typename ValueOf, typename Out>
std::uint32_t collectFlagged(std::uint32_t n, const std::uint32_t* flags,
                             const ValueOf& valueOf, Out* out) {
  std::uint32_t nKept = 0;
  for (std::uint32_t i = 0; i < n; ++i) {
    out[nKept] = valueOf(i);
    nKept += flags[i];
  }
  return nKept;
}

}  // namespace Acts::detail
