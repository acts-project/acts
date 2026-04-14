// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/CombinatoricIndices.hpp"

#include <ranges>
#include <set>

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
      BOOST_CHECK_EQUAL(indices.setSize(), n);
      BOOST_CHECK_EQUAL(indices.size(), binomial(n, 3ul));
    }
    if (n >= 2ul) {
      CombinatoricIndices<2> indices{n};
      BOOST_CHECK_EQUAL(indices.setSize(), n);
      BOOST_CHECK_EQUAL(indices.size(), binomial(n, 2ul));
    }
    CombinatoricIndices<1> indices{n};
    BOOST_CHECK_EQUAL(indices.setSize(), n);
    BOOST_CHECK_EQUAL(indices.size(), binomial(n, 1ul));
  }
}

BOOST_AUTO_TEST_CASE(CombinationDraw4D) {
  CombinatoricIndices<4> indexGenerator{10ul};
  BOOST_CHECK_EQUAL(indexGenerator.size(), binomial(10ul, 4ul));
  BOOST_CHECK_EQUAL(indexGenerator.setSize(), 10ul);

  CombinationSet<4> cachedCombos{};
  for (std::size_t comb = 0; comb < indexGenerator.size(); ++comb) {
    std::array<std::size_t, 4> combination = indexGenerator.draw(comb);
    std::ranges::sort(combination);
    // Check that all 4 indices are different
    for (std::size_t i = 1; i < 4; ++i) {
      for (std::size_t k = 0; k < i; ++k) {
        BOOST_CHECK_NE(combination[i], combination[k]);
      }
    }
    std::cout << "Iteration: " << comb << "-> drawn indices: " << combination
              << std::endl;
    BOOST_CHECK_EQUAL(cachedCombos.insert(combination).second, true);
  }
  BOOST_CHECK_EQUAL(cachedCombos.size(), indexGenerator.size());
}

BOOST_AUTO_TEST_SUITE_END()
