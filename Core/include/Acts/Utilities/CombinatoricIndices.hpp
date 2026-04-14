// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts {
/// Utilitty class to draw
template <std::size_t K>
class CombinatoricIndices {
 public:
  /// Declare other combinatoric indices classes as friends
  template <std::size_t L>
  friend class CombinatoricIndices;

  /// Constructor of the combinatoric indices
  /// @param The size of the set from which indices are drawn
  /// @note An exception is thrown if the size is less than the
  /// number of indices to draw
  explicit CombinatoricIndices(const std::size_t N);

  /// The number of possibilities to draw K different inidces from the set
  /// @returns The number of combinations
  std::size_t size() const;
  /// The size of the set from which the indices are drawn
  /// @returns The parameter N with which this class was constructed
  std::size_t setSize() const;

  /// Draws a new combination of indices and stores them in the array
  /// @param combination: The number of the combination in the sequence
  /// @returns An array where each is a unique number from [0 -N)
  std::array<std::size_t, K> draw(const std::size_t combination) const;

 private:
  std::size_t drawIndex(const std::size_t combination,
                        const std::size_t slot) const;

  static constexpr bool s_hasSubDraw = K > 1ul;
  using SubDrawer_t =
      std::conditional_t<s_hasSubDraw, CombinatoricIndices<K - 1>, std::size_t>;

  std::size_t m_N{};
  SubDrawer_t m_childDrawer{m_N - 1ul};

  std::vector<std::size_t> m_borders{};
};

}  // namespace Acts

#include "Acts/Utilities/CombinatoricIndices.ipp"
