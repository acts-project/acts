// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include <cmath>
#include <cstdint>

namespace Acts::detail {

/// Class that encapsulates a cantor pair, which represents an edge of a graph
/// By default ensures all edges are ordered, so the represented graph is
/// undirected: (a,b) and (b,a) are the same edge.
template <typename T>
class CantorEdge {
  T m_value;

 public:
  CantorEdge(T x, T y, bool sort = true) {
    if ((x > y) && sort) {
      std::swap(x, y);
    }
    m_value = y + ((x + y) * (x + y + 1)) / 2;
  }

  std::pair<T, T> inverse() const {
    auto f = [](T w) -> T { return (w * (w + 1)) / 2; };
    auto q = [](T w) -> T {
      return std::floor((std::sqrt(8 * w + 1) - 1) / 2);
    };

    auto y = m_value - f(q(m_value));
    auto x = q(m_value) - y;

    return {x, y};
  }

  T value() const { return m_value; }

  auto operator<=>(const CantorEdge<T>& other) const = default;
};

}  // namespace Acts::detail
