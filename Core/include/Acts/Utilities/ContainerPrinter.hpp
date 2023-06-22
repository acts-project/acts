// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <ostream>

namespace Acts {

template <typename iterator_t>
struct ContainerPrinter {
  iterator_t begin;
  iterator_t end;

  template <typename container_t>
  ContainerPrinter(const container_t &c, std::size_t max) {
    begin = c.cbegin();
    const auto n = std::min(c.size(), max);
    end = begin;
    std::advance(end, n);
    assert(std::distance(begin, end) == static_cast<std::ptrdiff_t>(n));
  }

  template <typename container_t>
  ContainerPrinter(const container_t &c) : ContainerPrinter(c, c.size()) {}

  ContainerPrinter(iterator_t a, iterator_t b) : begin(a), end(b) {}
};

template <typename iterator_t>
std::ostream &operator<<(std::ostream &os,
                         const ContainerPrinter<iterator_t> &p) {
  for (auto it = p.begin; it != p.end; ++it) {
    os << *it << " ";
  }
  return os;
}

template <class container_t>
ContainerPrinter(const container_t &, std::size_t)
    -> ContainerPrinter<typename std::decay_t<container_t>::const_iterator>;

template <class container_t>
ContainerPrinter(const container_t &)
    -> ContainerPrinter<typename std::decay_t<container_t>::const_iterator>;

}  // namespace Acts
