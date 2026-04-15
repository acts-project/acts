// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/CombinatoricIndices.hpp"

#include <ranges>
#include <set>

using namespace Acts::detail;
using namespace Acts;

struct ArraySorter {
  template <std::size_t K>
  bool operator()(const std::array<std::size_t, K>& a,
                  const std::array<std::size_t, K>& b) const {
    for (std::size_t i = 0; i < K; ++i) {
      if (a[i] != b[i]) {
        return a[i] < b[i];
      }
    }
    return false;
  }
};

template <std::size_t K>
std::ostream& operator<<(std::ostream& ostr,
                         const std::array<std::size_t, K>& x) {
  ostr << "[";
  for (std::size_t i = 0ul; i < x.size(); ++i) {
    ostr << x[i];
    if (i + 1ul != x.size()) {
      ostr << ", ";
    }
  }
  ostr << "]";
  return ostr;
}

template <std::size_t K>
using CombinationSet = std::set<std::array<std::size_t, K>, ArraySorter>;

BOOST_AUTO_TEST_SUITE(UtilitiesCloneablePtr)

BOOST_AUTO_TEST_CASE(SizeTest) {
  for (std::size_t n = 10ul; n >= 1ul; --n) {
    if (n >= 3ul) {
      CombinatoricIndices<3> indices{n};
      BOOST_CHECK_EQUAL(indices.intervalSize(), n);
      BOOST_CHECK_EQUAL(indices.size(), binomial(n, 3ul));
    }
    if (n >= 2ul) {
      CombinatoricIndices<2> indices{n};
      BOOST_CHECK_EQUAL(indices.intervalSize(), n);
      BOOST_CHECK_EQUAL(indices.size(), binomial(n, 2ul));
    }
    CombinatoricIndices<1> indices{n};
    BOOST_CHECK_EQUAL(indices.intervalSize(), n);
    BOOST_CHECK_EQUAL(indices.size(), binomial(n, 1ul));
  }
}

template <std::size_t K>
void checkComboDrawing(const std::size_t N) {
  if (N < K) {
    if constexpr (K > 1) {
      checkComboDrawing<K - 1>(N);
    }
    return;
  }
  const CombinatoricIndices<K> indexGenerator{N};
  BOOST_CHECK_EQUAL(indexGenerator.size(), binomial(N, K));
  BOOST_CHECK_EQUAL(indexGenerator.intervalSize(), N);

  CombinationSet<K> cachedCombos{};

  for (std::size_t combo = 0ul; combo < indexGenerator.size(); ++combo) {
    std::array<std::size_t, K> combination = indexGenerator.draw(combo);
    /// Check that all indices are less than N
    for (const std::size_t idx : combination) {
      BOOST_CHECK_LT(idx, N);
    }
    /// Check that all indices are unique
    for (std::size_t i = 1ul; i < combination.size(); ++i) {
      for (std::size_t k = 0ul; k < i; ++k) {
        BOOST_CHECK_LT(combination[k], combination[i]);
      }
    }
    /// Sort the indices as we are only inter
    std::ranges::sort(combination);
    std::cout << "Iteration: " << (combo) << "-> drawn indices: " << combination
              << std::endl;
    BOOST_CHECK_EQUAL(cachedCombos.insert(combination).second, true);
  }
  BOOST_CHECK_EQUAL(cachedCombos.size(), indexGenerator.size());
  if constexpr (K > 1) {
    checkComboDrawing<K - 1>(N);
  }
}

BOOST_AUTO_TEST_CASE(CombinationDraw) {
  checkComboDrawing<11ul>(11ul);
}

template <std::size_t K>
void checkCombinatoricIterator(const std::size_t N) {
  if (N < K) {
    return;
  }
  CombinatoricIndices<K> indexGenerator{N};
  BOOST_CHECK_EQUAL(indexGenerator.size(), binomial(N, K));
  BOOST_CHECK_EQUAL(indexGenerator.intervalSize(), N);

  auto begin = indexGenerator.begin();
  auto end = indexGenerator.end();
  BOOST_CHECK_EQUAL(begin != end, true);
  BOOST_CHECK_EQUAL((begin + indexGenerator.size()) == end, true);

  CombinationSet<K> cachedCombos{};
  for (const auto& combination : indexGenerator) {
    BOOST_CHECK_EQUAL(cachedCombos.insert(combination).second, true);

    for (const std::size_t idx : combination) {
      BOOST_CHECK_LT(idx, N);
    }
    /// Check that all indices are unique
    for (std::size_t i = 1ul; i < combination.size(); ++i) {
      for (std::size_t k = 0ul; k < i; ++k) {
        BOOST_CHECK_NE(combination[i], combination[k]);
      }
    }
  }
  if constexpr (K > 1) {
    checkCombinatoricIterator<K - 1>(N);
  }
}

BOOST_AUTO_TEST_CASE(CombinatoricIterator) {
  checkCombinatoricIterator<7>(16);
}

BOOST_AUTO_TEST_SUITE_END()
