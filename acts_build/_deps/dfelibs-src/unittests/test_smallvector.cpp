// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for dfe::SmallVector

#include <boost/test/unit_test.hpp>

#include <iterator>
#include <string>

#include "dfe/dfe_smallvector.hpp"

BOOST_AUTO_TEST_CASE(smallvector_int) {
  dfe::SmallVector<int, 4> sm;

  BOOST_TEST(sm.empty());
  BOOST_TEST(sm.size() == 0);
  BOOST_TEST(std::distance(sm.begin(), sm.end()) == 0);

  for (int i = 0; i < 1024; ++i) {
    BOOST_TEST(sm.emplace_back(i) == i);
    BOOST_TEST(sm.size() == (i + 1));
    BOOST_TEST(sm.size() == std::distance(sm.begin(), sm.end()));
    int j = 0;
    for (auto x : sm) {
      BOOST_TEST(sm[j] == j);
      BOOST_TEST(x == j);
      j += 1;
    }
    BOOST_TEST(j == sm.size());
  }
}

BOOST_AUTO_TEST_CASE(smallvector_int_emplace) {
  dfe::SmallVector<double, 6> at_front;
  dfe::SmallVector<double, 6> at_back;

  // insert same elements once at the back, once at the front
  for (int i = 0; i < 256; ++i) {
    double x = 1.23 * i;

    BOOST_TEST(*at_front.emplace(at_front.begin(), x) == x);
    BOOST_TEST(at_back.emplace_back(x) == x);
    for (int j = 0; j < (i + 1); ++j) {
      BOOST_TEST(at_front[j] == at_back[i - j]);
    }
  }
}

class Simple {
public:
  Simple() = default;
  Simple(int i) : m_num(2 * i), m_str("abc") {
    for (; 0 < i; --i) {
      m_str.append("xy");
    }
  }

  int num() const { return m_num; }
  const std::string str() const { return m_str; }

private:
  int m_num = 42;
  std::string m_str;
};

BOOST_AUTO_TEST_CASE(smallvector_simpleclass) {
  dfe::SmallVector<Simple, 8> sm;

  BOOST_TEST(sm.empty());
  BOOST_TEST(sm.size() == 0);
  BOOST_TEST(std::distance(sm.begin(), sm.end()) == 0);

  for (int i = 0; i < 128; ++i) {
    BOOST_TEST(sm.emplace_back(i).num() == (2 * i));
    BOOST_TEST(sm.size() == (i + 1));
    BOOST_TEST(sm.size() == std::distance(sm.begin(), sm.end()));
    int j = 0;
    for (auto x : sm) {
      Simple cmp(j);
      BOOST_TEST(sm[j].num() == cmp.num());
      BOOST_TEST(sm[j].str() == cmp.str());
      BOOST_TEST(x.num() == cmp.num());
      BOOST_TEST(x.str() == cmp.str());
      j += 1;
    }
    BOOST_TEST(j == sm.size());
  }
}
