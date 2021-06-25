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
    BOOST_CHECK_EQUAL(pv.get(), nullptr);

    BOOST_CHECK(IsPolymorphicValue<decltype(pv)>::value);
    BOOST_CHECK(!IsPolymorphicValue<A>::value);
  }

  {
    PolymorphicValue<Base> pv;
    BOOST_CHECK(!pv);
    pv = A{6};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_NE(pv.get(), nullptr);
    BOOST_CHECK_EQUAL(pv->value(), 6);
    BOOST_CHECK_EQUAL((*pv).value(), 6);
  }

  {
    PolymorphicValue<Base> pv{A{8}};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_NE(pv.get(), nullptr);
    BOOST_CHECK_EQUAL(pv->value(), 8);

    pv = A{9};
    BOOST_CHECK(!!pv);
    BOOST_CHECK_NE(pv.get(), nullptr);
    BOOST_CHECK_EQUAL(pv->value(), 9);
  }
}

BOOST_AUTO_TEST_CASE(TestMakePolymorphicValue) {
  {
    PolymorphicValue<Base> pv{makePolymorphicValue<Base, A>(19)};
    BOOST_CHECK(pv);
    BOOST_CHECK_EQUAL(pv->value(), 19);
  }
  {
    PolymorphicValue<Base> pv;
    BOOST_CHECK(!pv);
    pv = makePolymorphicValue<A>(37);
    BOOST_CHECK(pv);
    BOOST_CHECK_EQUAL(pv->value(), 37);
  }
}

struct Destruct : public Base {
  bool* m_destroyed;
  int m_value;
  Destruct(bool* d, int v = 0) : m_destroyed{d}, m_value{v} {}

  int value() override { return m_value; }
  ~Destruct() { (*m_destroyed) = true; }
};

BOOST_AUTO_TEST_CASE(TestConstructFromPointer) {
  {
    bool destroyed = false;
    Destruct* p = new Destruct(&destroyed, 62);

    PolymorphicValue<Destruct> pv{p};
    BOOST_CHECK(pv);
    BOOST_CHECK(!destroyed);
    BOOST_CHECK_EQUAL(pv->value(), 62);
    BOOST_CHECK_EQUAL(pv.get(), p);

    auto pv2 = pv;  // copy
    BOOST_CHECK_NE(pv2.get(), pv.get());

    pv.reset();
    BOOST_CHECK(!pv);
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(pv.get(), nullptr);
  }

  {
    bool destroyed = false;
    // this is pointer to base, we can't copy this
    Base* p = new Destruct(&destroyed, 62);

    PolymorphicValue<Base> pv{p};
    BOOST_CHECK(pv);
    BOOST_CHECK(!destroyed);
    BOOST_CHECK_EQUAL(pv->value(), 62);
    BOOST_CHECK_EQUAL(pv.get(), p);

    // copy is not possible, when PV was constructed, we only got Base* so we
    // couldn't figure out how copying works
    BOOST_CHECK_THROW(auto pv2 = pv, std::runtime_error);

    pv.reset();
    BOOST_CHECK(!pv);
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(pv.get(), nullptr);
  }
}

BOOST_AUTO_TEST_CASE(TestDestruction) {
  {
    bool destroyed = false;

    std::optional<PolymorphicValue<Base>> pvOpt =
        makePolymorphicValue<Destruct>(&destroyed);

    BOOST_CHECK(!destroyed);
    BOOST_CHECK_EQUAL((*pvOpt)->value(), 0);
    pvOpt = std::nullopt;
    BOOST_CHECK(!!destroyed);
  }

  {
    bool destroyed = false;

    PolymorphicValue<Base> pv{makePolymorphicValue<Destruct>(&destroyed)};

    BOOST_CHECK(!destroyed);
    BOOST_CHECK_EQUAL(pv->value(), 0);
    pv.reset();
    BOOST_CHECK(!!destroyed);
    BOOST_CHECK(!pv);
  }
}

BOOST_AUTO_TEST_CASE(TestRelease) {
  {
    bool destroyed = false;

    std::optional<PolymorphicValue<Base>> pvOpt =
        makePolymorphicValue<Destruct>(&destroyed);

    BOOST_CHECK(!destroyed);
    BOOST_CHECK_EQUAL((*pvOpt)->value(), 0);

    BOOST_CHECK((*pvOpt));
    BOOST_CHECK_NE((*pvOpt).get(), nullptr);
    Base* p = (*pvOpt).release();
    BOOST_CHECK(!(*pvOpt));
    BOOST_CHECK_EQUAL((*pvOpt).get(), nullptr);

    pvOpt = std::nullopt;
    BOOST_CHECK(!destroyed);

    delete p;
    BOOST_CHECK(destroyed);
  }
}

struct Base1 {
  virtual ~Base1() = default;
  virtual int value() = 0;
};

struct E : public Base1 {
  E(int v) : m_value(v) {}
  int value() override { return m_value; }
  int m_value;
};

struct F : public Base1 {
  F(int v) : m_value(v) {}
  int value() override { return m_value; }
  int m_value;
};

BOOST_AUTO_TEST_CASE(TestAssignmentSameBase) {
  {  // COPY
    PolymorphicValue<Base1> pv1{E{8}};
    PolymorphicValue<Base1> pv2{F{2}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<F*>(pv1.get()), nullptr);

    BOOST_CHECK_EQUAL(pv2->value(), 2);
    BOOST_CHECK_NE(dynamic_cast<F*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());

    Base1* pre = pv2.get();

    pv1 = pv2;

    BOOST_CHECK(pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 2);
    BOOST_CHECK_NE(dynamic_cast<F*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_EQUAL(pv1->value(), 2);
    BOOST_CHECK_NE(dynamic_cast<F*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv1.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);

    BOOST_CHECK_NE(pv1.get(), pre);
  }

  {  // MOVE
    PolymorphicValue<Base1> pv1{E{8}};
    PolymorphicValue<Base1> pv2{F{2}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<F*>(pv1.get()), nullptr);

    BOOST_CHECK_EQUAL(pv2->value(), 2);
    BOOST_CHECK_NE(dynamic_cast<F*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());

    Base1* pre = pv2.get();

    pv1 = std::move(pv2);

    BOOST_CHECK(pv1);
    BOOST_CHECK(!pv2);

    BOOST_CHECK_EQUAL(pv1->value(), 2);
    BOOST_CHECK_NE(dynamic_cast<F*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv1.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_EQUAL(pv2.get(), nullptr);

    BOOST_CHECK_EQUAL(pv1.get(), pre);
  }
}

BOOST_AUTO_TEST_CASE(TestConstructionSameBase) {
  {  // COPY
    PolymorphicValue<Base1> pv1{E{8}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);

    PolymorphicValue<Base1> pv2{pv1};

    BOOST_CHECK(pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);
  }

  {  // MOVE
    PolymorphicValue<Base1> pv1{E{8}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);

    PolymorphicValue<Base1> pv2{std::move(pv1)};

    BOOST_CHECK(!pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_EQUAL(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);
  }
}

struct Base2 : public Base1 {};
struct G : public Base2 {
  G(int v) : m_value(v) {}
  int value() override { return m_value; }
  int m_value;
};

BOOST_AUTO_TEST_CASE(TestAssignmentDifferentBase) {
  {  // COPY
    PolymorphicValue<Base1> pv1{E{8}};
    PolymorphicValue<Base2> pv2{G{3}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<G*>(pv1.get()), nullptr);

    BOOST_CHECK_EQUAL(pv2->value(), 3);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());

    Base1* pre = pv2.get();
    pv1 = pv2;

    BOOST_CHECK(pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 3);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_EQUAL(pv1->value(), 3);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv1.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);

    BOOST_CHECK_NE(pv1.get(), pre);
  }

  {  // MOVE
    PolymorphicValue<Base1> pv1{E{8}};
    PolymorphicValue<Base2> pv2{G{3}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<E*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<G*>(pv1.get()), nullptr);

    BOOST_CHECK_EQUAL(pv2->value(), 3);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv2.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());

    Base1* pre = pv2.get();
    pv1 = std::move(pv2);

    BOOST_CHECK(pv1);
    BOOST_CHECK(!pv2);

    BOOST_CHECK_EQUAL(pv1->value(), 3);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv1.get()), nullptr);
    BOOST_CHECK_EQUAL(dynamic_cast<E*>(pv1.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_EQUAL(pv2.get(), nullptr);

    BOOST_CHECK_EQUAL(pv1.get(), pre);
  }
}

BOOST_AUTO_TEST_CASE(TestConstructionDifferentBase) {
  {  // COPY
    BOOST_CHECK(std::is_copy_constructible_v<G>);
    PolymorphicValue<Base2> pv1{G{8}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv1.get()), nullptr);

    PolymorphicValue<Base1> pv2{pv1};

    BOOST_CHECK(pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_NE(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);
  }

  {  // MOVE
    PolymorphicValue<Base2> pv1{G{8}};

    BOOST_CHECK_EQUAL(pv1->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv1.get()), nullptr);

    PolymorphicValue<Base1> pv2{std::move(pv1)};

    BOOST_CHECK(!pv1);
    BOOST_CHECK(pv2);

    BOOST_CHECK_EQUAL(pv2->value(), 8);
    BOOST_CHECK_NE(dynamic_cast<G*>(pv2.get()), nullptr);

    BOOST_CHECK_NE(pv1.get(), pv2.get());
    BOOST_CHECK_EQUAL(pv1.get(), nullptr);
    BOOST_CHECK_NE(pv2.get(), nullptr);
  }
}

struct Destruct2 : public Base2 {
  bool* m_destroyed;
  Destruct2(bool* d) : m_destroyed{d} {}

  int value() override { return 0; }
  ~Destruct2() { (*m_destroyed) = true; }
};

BOOST_AUTO_TEST_CASE(TestDestroyDelegate) {
  bool destroyed = false;
  std::optional<PolymorphicValue<Base2>> pvOpt2{
      makePolymorphicValue<Destruct2>(&destroyed)};

  BOOST_CHECK(!destroyed);
  BOOST_CHECK_EQUAL((*pvOpt2)->value(), 0);

  std::optional<PolymorphicValue<Base1>> pvOpt1{std::move(*pvOpt2)};
  BOOST_CHECK(!destroyed);

  pvOpt2 = std::nullopt;
  BOOST_CHECK(!destroyed);

  BOOST_CHECK_EQUAL((*pvOpt1)->value(), 0);

  pvOpt1 = std::nullopt;
  BOOST_CHECK(destroyed);
}

BOOST_AUTO_TEST_CASE(TestReleaseDelegate) {
  bool destroyed = false;
  std::optional<PolymorphicValue<Base2>> pvOpt2{
      makePolymorphicValue<Destruct2>(&destroyed)};

  BOOST_CHECK(!destroyed);
  BOOST_CHECK_EQUAL((*pvOpt2)->value(), 0);

  std::optional<PolymorphicValue<Base1>> pvOpt1{std::move(*pvOpt2)};
  BOOST_CHECK(!destroyed);

  pvOpt2 = std::nullopt;
  BOOST_CHECK(!destroyed);

  BOOST_CHECK_EQUAL((*pvOpt1)->value(), 0);

  BOOST_CHECK((*pvOpt1));
  BOOST_CHECK_NE((*pvOpt1).get(), nullptr);
  Base1* p = (*pvOpt1).release();
  BOOST_CHECK(!(*pvOpt1));
  BOOST_CHECK_EQUAL((*pvOpt1).get(), nullptr);

  pvOpt1 = std::nullopt;
  BOOST_CHECK(!destroyed);

  delete p;
  BOOST_CHECK(destroyed);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace Acts
