// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <cstdint>

#include <boost/safe_numerics/safe_integer.hpp>

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

    // Use safe integers here to ensure that no overflow happens
    using namespace boost::safe_numerics;

    safe<T> sx(x);
    safe<T> sy(y);

    auto res = sy + ((sx + sy) * (sx + sy + 1)) / 2;

    m_value = res;
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
