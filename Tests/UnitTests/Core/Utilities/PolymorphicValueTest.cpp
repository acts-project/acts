// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/PolymorphicValue.hpp"

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Utilities)

struct Base {
  virtual ~Base() = default;
  virtual int value() = 0;
};

struct A : public Base {
  A(int v) : m_value{v} {}

  int m_value;
  int value() override { return m_value; }
};

struct D {};

BOOST_AUTO_TEST_CASE(TestConstruction) {
  {
    PolymorphicValue<Base> pv;
    BOOST_CHECK(!pv);
    // BOOST_CHECK_THROW(PolymorphicValue<Base>{D{}}, std::invalid_argument);
    // BOOST_CHECK_THROW(pv = D{}, std::invalid_argument);
  }

  {
    PolymorphicValue<Base> pv;
    BOOST_CHECK(!pv);
    pv = A{6};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_EQUAL(pv->value(), 6);
  }

  {
    PolymorphicValue<Base> pv{A{8}};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_EQUAL(pv->value(), 8);

    pv = A{9};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_EQUAL(pv->value(), 9);
  }
}

struct Destruct : public Base {
  bool* m_destroyed;
  Destruct(bool* d) : m_destroyed{d} {}

  int value() override { return 0; }
  ~Destruct() { (*m_destroyed) = true; }
};

BOOST_AUTO_TEST_CASE(TestDestruction) {
  {
    bool destroyed = false;

    std::optional<PolymorphicValue<Base>> pvOpt =
        PolymorphicValue<Base>{std::in_place_type_t<Destruct>(), &destroyed};

    BOOST_CHECK(!destroyed);
    pvOpt = std::nullopt;
    BOOST_CHECK(!!destroyed);
  }

  {
    bool destroyed = false;

    PolymorphicValue<Base> pv{std::in_place_type_t<Destruct>(), &destroyed};

    BOOST_CHECK(!destroyed);
    pv.reset();
    BOOST_CHECK(!!destroyed);
    BOOST_CHECK(!pv);
  }
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
