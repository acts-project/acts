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

template <std::size_t K>
bool goodArray(const std::array<std::size_t, K>& array, const std::size_t N) {
  bool good{true};
  for (std::size_t k = 0; k < array.size(); ++k) {
    BOOST_CHECK_LT(array[k], N);
    good &= (array[k] < N);
  }
  for (std::size_t j = 1; j < array.size(); ++j) {
    for (std::size_t i = 0; i < j; ++i) {
      BOOST_CHECK_NE(array[j], array[i]);
      good &= (array[i] != array[j]);
    }
  }
  return good;
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
    std::cout << "N=" << N << ", K=" << K << " --- Iteration: " << (combo)
              << "-> drawn indices: " << combination << std::endl;
    BOOST_CHECK_EQUAL(goodArray(combination, N), true);
    /// Check that all indices are sorted
    for (std::size_t i = 1ul; i < combination.size(); ++i) {
      for (std::size_t k = 0ul; k < i; ++k) {
        BOOST_CHECK_LT(combination[k], combination[i]);
      }
    }
    BOOST_CHECK_EQUAL(cachedCombos.insert(combination).second, true);
  }
  BOOST_CHECK_EQUAL(cachedCombos.size(), indexGenerator.size());
  if constexpr (K > 1) {
    checkComboDrawing<K - 1>(N);
  }
}

BOOST_AUTO_TEST_CASE(CombinationDraw) {
  return;
  checkComboDrawing<11ul>(11ul);
  checkComboDrawing<10ul>(10ul);
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
    BOOST_CHECK_EQUAL(goodArray(combination, N), true);
    BOOST_CHECK_EQUAL(cachedCombos.insert(combination).second, true);
  }
  if constexpr (K > 1) {
    checkCombinatoricIterator<K - 1>(N);
  }
}

BOOST_AUTO_TEST_CASE(CombinatoricIterator) {
  return;
  checkCombinatoricIterator<7>(16);
}

BOOST_AUTO_TEST_CASE(RemapIndices4) {
  constexpr std::size_t N = 10;
  CombinatoricIndices<4> indexGenerator{N};
  CombinationSet<4> cachedCombos{}, allCombos{};

  for (std::size_t combo = 0; combo < indexGenerator.size(); ++combo) {
    const auto indices = indexGenerator.draw(combo);
    allCombos.insert(indices);
    auto indexCopy = indices;
    indexCopy[3] = N - indices[3] + indices[1];
    indexCopy[2] = N - indices[2] + indices[1];

    std::cout << "RemapIndices - Iteration: " << (combo + 1)
              << " -> drawn: " << indices << " --> reshuffled: " << indexCopy
              << std::endl;
    const bool goodone = goodArray(indexCopy, N);
    BOOST_CHECK_EQUAL(goodone, true);
    std::ranges::sort(indexCopy);
    if (goodone) {
      BOOST_CHECK_EQUAL(cachedCombos.insert(indexCopy).second, true);
    }
  }
  BOOST_CHECK_EQUAL(allCombos.size(), cachedCombos.size());
  for (const auto& combo : allCombos) {
    if (!cachedCombos.count(combo)) {
      std::cout << "ReMapIndices -- Combination did not made it " << combo
                << std::endl;
    }
  }
}
BOOST_AUTO_TEST_CASE(RemapIndices3) {
  constexpr std::size_t N = 10;
  CombinatoricIndices<3> indexGenerator{N};
  CombinationSet<3> cachedCombos{}, allCombos{};

  for (std::size_t combo = 0; combo < indexGenerator.size(); ++combo) {
    const auto indices = indexGenerator.draw(combo);
    allCombos.insert(indices);
    auto indexCopy = indices;
    indexCopy[2] = N - indices[2] + indices[1];

    std::cout << "RemapIndices - Iteration: " << (combo + 1)
              << " -> drawn: " << indices << " --> reshuffled: " << indexCopy
              << std::endl;
    const bool goodone = goodArray(indexCopy, N);
    BOOST_CHECK_EQUAL(goodone, true);
    std::ranges::sort(indexCopy);
    if (goodone) {
      BOOST_CHECK_EQUAL(cachedCombos.insert(indexCopy).second, true);
    }
  }
  BOOST_CHECK_EQUAL(allCombos.size(), cachedCombos.size());
  for (const auto& combo : allCombos) {
    if (!cachedCombos.count(combo)) {
      std::cout << "ReMapIndices -- Combination did not made it " << combo
                << std::endl;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
