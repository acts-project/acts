// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/CloneablePtr.hpp"

#include <memory>
#include <string>
#include <utility>

using namespace Acts;

namespace {

struct Base {
  virtual ~Base() = default;
  virtual int value() const = 0;
  virtual std::unique_ptr<Base> clone() const = 0;
};

struct Derived : Base {
  int m_val;
  explicit Derived(int v) : m_val(v) {}
  int value() const override { return m_val; }
  std::unique_ptr<Base> clone() const override {
    return std::make_unique<Derived>(m_val);
  }
};

}  // namespace

BOOST_AUTO_TEST_SUITE(UtilitiesCloneablePtr)

BOOST_AUTO_TEST_CASE(DefaultConstruction) {
  CloneablePtr<int> p;
  BOOST_CHECK(!p);
  BOOST_CHECK(p.get() == nullptr);
}

BOOST_AUTO_TEST_CASE(ConstructFromUniquePtr) {
  CloneablePtr<int> p(std::make_unique<int>(42));
  BOOST_CHECK(p);
  BOOST_CHECK_EQUAL(*p, 42);
}

BOOST_AUTO_TEST_CASE(ConstructFromRawPtr) {
  CloneablePtr<int> p(new int(7));
  BOOST_CHECK(p);
  BOOST_CHECK_EQUAL(*p, 7);
}

BOOST_AUTO_TEST_CASE(CopyConstruction) {
  CloneablePtr<int> a(std::make_unique<int>(10));
  CloneablePtr<int> b(a);
  BOOST_CHECK(a);
  BOOST_CHECK(b);
  BOOST_CHECK_EQUAL(*a, 10);
  BOOST_CHECK_EQUAL(*b, 10);
  BOOST_CHECK(a.get() != b.get());
}

BOOST_AUTO_TEST_CASE(CopyAssignment) {
  CloneablePtr<int> a(std::make_unique<int>(20));
  CloneablePtr<int> b;
  b = a;
  BOOST_CHECK(b);
  BOOST_CHECK_EQUAL(*b, 20);
  BOOST_CHECK(a.get() != b.get());
}

BOOST_AUTO_TEST_CASE(MoveConstruction) {
  CloneablePtr<int> a(std::make_unique<int>(30));
  int* raw = a.get();
  CloneablePtr<int> b(std::move(a));
  BOOST_CHECK(b);
  BOOST_CHECK_EQUAL(b.get(), raw);
  BOOST_CHECK(!a);  // NOLINT(bugprone-use-after-move)
}

BOOST_AUTO_TEST_CASE(MoveAssignment) {
  CloneablePtr<int> a(std::make_unique<int>(40));
  int* raw = a.get();
  CloneablePtr<int> b;
  b = std::move(a);
  BOOST_CHECK(b);
  BOOST_CHECK_EQUAL(b.get(), raw);
  BOOST_CHECK(!a);  // NOLINT(bugprone-use-after-move)
}

BOOST_AUTO_TEST_CASE(NullCopy) {
  CloneablePtr<int> a;
  CloneablePtr<int> b(a);
  BOOST_CHECK(!a);
  BOOST_CHECK(!b);
}

BOOST_AUTO_TEST_CASE(Release) {
  CloneablePtr<int> p(std::make_unique<int>(50));
  int* raw = p.release();
  BOOST_CHECK(!p);
  BOOST_CHECK_EQUAL(*raw, 50);
  delete raw;
}

BOOST_AUTO_TEST_CASE(Reset) {
  CloneablePtr<int> p(std::make_unique<int>(60));
  p.reset();
  BOOST_CHECK(!p);

  p.reset(new int(70));
  BOOST_CHECK(p);
  BOOST_CHECK_EQUAL(*p, 70);
}

BOOST_AUTO_TEST_CASE(ArrowOperator) {
  CloneablePtr<std::string> p(std::make_unique<std::string>("hello"));
  BOOST_CHECK_EQUAL(p->size(), 5u);
}

BOOST_AUTO_TEST_CASE(CustomCloner) {
  auto cloner = [](const Base& b) { return b.clone(); };

  CloneablePtr<Base> a(std::make_unique<Derived>(99), std::move(cloner));
  BOOST_CHECK(a);
  BOOST_CHECK_EQUAL(a->value(), 99);

  CloneablePtr<Base> b(a);
  BOOST_CHECK(b);
  BOOST_CHECK_EQUAL(b->value(), 99);
  BOOST_CHECK(a.get() != b.get());
}

BOOST_AUTO_TEST_CASE(NullptrComparison) {
  CloneablePtr<int> null;
  CloneablePtr<int> nonNull(std::make_unique<int>(1));

  BOOST_CHECK(null == nullptr);
  BOOST_CHECK(nullptr == null);
  BOOST_CHECK(!(null != nullptr));

  BOOST_CHECK(nonNull != nullptr);
  BOOST_CHECK(nullptr != nonNull);
  BOOST_CHECK(!(nonNull == nullptr));
}

BOOST_AUTO_TEST_CASE(SelfAssignment) {
  CloneablePtr<int> p(std::make_unique<int>(5));
  auto* addr = p.get();
  auto& ref = p;
  p = ref;
  BOOST_CHECK(p);
  BOOST_CHECK_EQUAL(*p, 5);
  BOOST_CHECK_EQUAL(p.get(), addr);
}

BOOST_AUTO_TEST_SUITE_END()
