// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/BFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace tt = boost::test_tools;

namespace Acts {

namespace Test {

// Create a test context
MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_CASE(TypeErasedCacheType) {
  bool constructor_called = false;
  bool destructor_called = false;

  struct MyCache {
    MyCache(int value, bool* ctor, bool* dtor) : m_value{value}, m_dtor{dtor} {
      (*ctor) = true;
    }
    ~MyCache() { (*m_dtor) = true; }
    int m_value;
    bool* m_dtor;
  };

  BOOST_CHECK(!constructor_called);
  BOOST_CHECK(!destructor_called);

  {
    BFieldProvider::Cache cache{BFieldProvider::Cache::make<MyCache>(
        42, &constructor_called, &destructor_called)};
    BOOST_CHECK(constructor_called);
    BOOST_CHECK(!destructor_called);

    MyCache& v = cache.get<MyCache>();
    BOOST_CHECK_EQUAL(v.m_value, 42);
    v.m_value = 65;

    MyCache& v2 = cache.get<MyCache>();
    BOOST_CHECK_EQUAL(v2.m_value, 65);
  }

  BOOST_CHECK(constructor_called);
  BOOST_CHECK(destructor_called);
}

BOOST_AUTO_TEST_CASE(CacheNonCopyable) {
  struct MyCache {
    int m_value{0};

    MyCache() = default;
    MyCache(int value) : m_value(value) {}

    MyCache(const MyCache&) = delete;
    MyCache& operator=(const MyCache&) = delete;

    MyCache& operator=(MyCache&&) = default;
    MyCache(MyCache&&) = default;
  };

  auto cache = BFieldProvider::Cache::make<MyCache>(42);
  // MyCache c{42};
  // std::any a;
  // a = std::move(c);
  // a.emplace<MyCache>(42);
  // a = std::make_any<MyCache>(42);
}

}  // namespace Test

}  // namespace Acts
