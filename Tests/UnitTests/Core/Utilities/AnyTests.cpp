// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Any.hpp"

#include <utility>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(AnyTests)

BOOST_AUTO_TEST_CASE(AnyConstructPrimitive) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL(a.as<int>(), v);
    BOOST_CHECK_NE(a.as<int>(), v + 1);
  }

  {
    // type that is large
    Any a;
    BOOST_CHECK(!a);

    std::array<int, 2> v{1, 2};
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<decltype(v)>().begin(),
                                  a.as<decltype(v)>().end(), v.begin(),
                                  v.end());
  }

  {
    // type that is large
    Any a;
    BOOST_CHECK(!a);

    std::array<unsigned long, 5> v{1, 2, 3, 4, 5};
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<decltype(v)>().begin(),
                                  a.as<decltype(v)>().end(), v.begin(),
                                  v.end());
  }
}

BOOST_AUTO_TEST_CASE(AnyConstructCustom) {
  struct A {
    int value;
    A() { value = 76; }
  };

  Any a;
  BOOST_CHECK(!a);
  a = Any{A{}};

  BOOST_CHECK(!!a);

  BOOST_CHECK_EQUAL(a.as<A>().value, 76);
}

BOOST_AUTO_TEST_CASE(AnyConstructCustomInPlace) {
  struct A {
    int value;
    A(int v) { value = v; }
  };

  Any a{std::in_place_type<A>, 42};
  BOOST_CHECK(!!a);
  BOOST_CHECK_EQUAL(a.as<A>().value, 42);
}

BOOST_AUTO_TEST_CASE(AnyMove) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    Any b = std::move(a);
    BOOST_CHECK(!!b);
    BOOST_CHECK_EQUAL(b.as<int>(), 5);

    Any c;
    c = std::move(b);
    BOOST_CHECK(!!c);
    BOOST_CHECK_EQUAL(c.as<int>(), 5);
  }
}

BOOST_AUTO_TEST_CASE(AnyCopy) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    Any b = a;
    BOOST_CHECK(!!b);
    BOOST_CHECK_EQUAL(b.as<int>(), 5);

    Any c;
    c = a;
    BOOST_CHECK(!!c);
    BOOST_CHECK_EQUAL(c.as<int>(), 5);
  }
}

struct D {
  bool* destroyed;
  D(bool* d) : destroyed{d} {}
  ~D() { *destroyed = true; }
};

struct D2 {
  bool* destroyed{nullptr};
  std::array<char, 512> blob{};

  D2(bool* d) : destroyed{d} {}

  ~D2() { *destroyed = true; }
};

BOOST_AUTO_TEST_CASE(AnyDestroy) {
  {  // small type
    bool destroyed = false;
    D d{&destroyed};
    BOOST_CHECK(!destroyed);
    {
      Any a{std::move(d)};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }

  {  // large type
    bool destroyed = false;
    D2 d{&destroyed};
    BOOST_CHECK(!destroyed);
    {
      Any a{std::move(d)};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
}

BOOST_AUTO_TEST_CASE(AnyDestroyInPlace) {
  {  // small type
    bool destroyed = false;
    BOOST_CHECK(!destroyed);
    {
      Any a{std::in_place_type<D>, &destroyed};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }

  {  // large type
    bool destroyed = false;
    BOOST_CHECK(!destroyed);
    {
      Any a{std::in_place_type<D2>, &destroyed};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
}

struct D3 {
  size_t* destroyed{nullptr};
  std::array<char, 512> blob{};

  D3(size_t* d) : destroyed{d} {}

  ~D3() { (*destroyed)++; }
};

BOOST_AUTO_TEST_CASE(LeakCheck) {
  size_t destroyed = 0;
  for (size_t i = 0; i < 10000; i++) {
    {
      BOOST_CHECK_EQUAL(destroyed, i);
      Any a;
      BOOST_CHECK_EQUAL(destroyed, i);
      a = Any{std::in_place_type<D3>, &destroyed};
      BOOST_CHECK_EQUAL(destroyed, i);
    }
    BOOST_CHECK_EQUAL(destroyed, i + 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()
