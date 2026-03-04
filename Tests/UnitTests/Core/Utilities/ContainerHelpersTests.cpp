// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/ContainerIterator.hpp"
#include "Acts/Utilities/detail/ContainerRange.hpp"
#include "Acts/Utilities/detail/ContainerSubset.hpp"

#include <set>

using namespace Acts;
using namespace Acts::detail;

namespace ActsTests {

class Container {
 public:
  Container() = default;

  Container(std::initializer_list<int> init) : m_data(init) {}

  int& operator[](std::size_t index) { return m_data[index]; }

  int operator[](std::size_t index) const { return m_data[index]; }

  int& at(std::size_t index) {
    if (index >= m_data.size()) {
      throw std::out_of_range("Index out of bounds");
    }
    return m_data.at(index);
  }

  int at(std::size_t index) const {
    if (index >= m_data.size()) {
      throw std::out_of_range("Index out of bounds");
    }
    return m_data.at(index);
  }

  template <bool ReadOnly>
  using Iterator = ContainerIterator<Container, int, std::size_t, ReadOnly>;

  Iterator<false> begin() { return {*this, 0}; }

  Iterator<false> end() { return {*this, m_data.size()}; }

  Iterator<true> begin() const { return {*this, 0}; }

  Iterator<true> end() const { return {*this, m_data.size()}; }

  template <bool ReadOnly>
  class Range : public ContainerRange<Range<ReadOnly>, Range<true>, Container,
                                      std::size_t, ReadOnly> {
   public:
    using Base = ContainerRange<Range<ReadOnly>, Range<true>, Container,
                                std::size_t, ReadOnly>;

    using Base::Base;
  };

  Range<false> range(std::size_t offset, std::size_t count) {
    return {*this, {offset, count}};
  }

  Range<true> range(std::size_t offset, std::size_t count) const {
    return {*this, {offset, count}};
  }

  template <bool ReadOnly>
  class Subset
      : public ContainerSubset<Subset<ReadOnly>, Subset<true>, Container, int,
                               std::span<const std::size_t>, ReadOnly> {
   public:
    using Base = ContainerSubset<Subset<ReadOnly>, Subset<true>, Container, int,
                                 std::span<const std::size_t>, ReadOnly>;

    using Base::Base;
  };

  Subset<false> subset(std::span<const std::size_t> indices) {
    return {*this, indices};
  }

  Subset<true> subset(std::span<const std::size_t> indices) const {
    return {*this, indices};
  }

  template <bool ReadOnly>
  class OwningSubset
      : public ContainerSubset<OwningSubset<ReadOnly>, OwningSubset<true>,
                               Container, int, std::vector<std::size_t>,
                               ReadOnly> {
   public:
    using Base =
        ContainerSubset<OwningSubset<ReadOnly>, OwningSubset<true>, Container,
                        int, std::vector<std::size_t>, ReadOnly>;

    using Base::Base;
  };

  OwningSubset<false> owningSubset(std::vector<std::size_t> indices) {
    return {*this, std::move(indices)};
  }

  OwningSubset<true> owningSubset(std::vector<std::size_t> indices) const {
    return {*this, std::move(indices)};
  }

  template <bool ReadOnly>
  class OwningOrderedSubset
      : public ContainerSubset<OwningOrderedSubset<ReadOnly>,
                               OwningOrderedSubset<true>, Container, int,
                               std::set<std::size_t>, ReadOnly> {
   public:
    using Base = ContainerSubset<OwningOrderedSubset<ReadOnly>,
                                 OwningOrderedSubset<true>, Container, int,
                                 std::set<std::size_t>, ReadOnly>;

    using Base::Base;
  };

  OwningOrderedSubset<false> owningOrderedSubset(
      std::set<std::size_t> indices) {
    return {*this, std::move(indices)};
  }

  OwningOrderedSubset<true> owningOrderedSubset(
      std::set<std::size_t> indices) const {
    return {*this, std::move(indices)};
  }

 private:
  std::vector<int> m_data;
};

BOOST_AUTO_TEST_SUITE(ContainerHelpersTests)

BOOST_AUTO_TEST_CASE(All) {
  Container mutableContainer{10, 20, 30, 40, 50};
  const Container& constContainer = mutableContainer;

  {
    std::vector<int> expected{10, 20, 30, 40, 50};
    std::vector<int> actual;
    for (const int i : mutableContainer) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<int> expected{10, 20, 30, 40, 50};
    std::vector<int> actual;
    for (const int i : constContainer) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<int> expected{20, 30, 40};
    std::vector<int> actual;
    for (const int i : mutableContainer.range(1, 3)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<int> expected{20, 30, 40};
    std::vector<int> actual;
    for (const int i : constContainer.range(1, 3)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<std::size_t> indices{0, 2, 4};
    std::vector<int> expected{10, 30, 50};
    std::vector<int> actual;
    for (const int i : constContainer.subset(indices)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<std::size_t> indices{0, 2, 4};
    std::vector<int> expected{10, 30, 50};
    std::vector<int> actual;
    for (const int i : mutableContainer.owningSubset(indices)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::vector<std::size_t> indices{0, 2, 4};
    std::vector<int> expected{10, 30, 50};
    std::vector<int> actual;
    for (const int i : constContainer.owningSubset(indices)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::set<std::size_t> indices{0, 2, 4};
    std::vector<int> expected{10, 30, 50};
    std::vector<int> actual;
    for (const int i : mutableContainer.owningOrderedSubset(indices)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }

  {
    std::set<std::size_t> indices{0, 2, 4};
    std::vector<int> expected{10, 30, 50};
    std::vector<int> actual;
    for (const int i : constContainer.owningOrderedSubset(indices)) {
      actual.push_back(i);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(expected.begin(), expected.end(),
                                  actual.begin(), actual.end());
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
