// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/TypeTraits.hpp"

#include <type_traits>

using namespace Acts::Concepts;

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

// generate traits for methods named foo and bar
METHOD_TRAIT(foo_method_t, foo);
METHOD_TRAIT(bar_method_t, bar);

struct E {
  int bar(const double& /*unused*/) { return 5; }
};

struct E2 {
  int bar(const double& /*unused*/) const { return 5; }
};

class E3 {
  int bar(const double& /*unused*/) { return 5; }
};

BOOST_AUTO_TEST_CASE(TypeTraitsMethods) {
  // E does not have a method bar without arguments
  static_assert(!has_method<E, int, bar_method_t>, "failed");
  // E does have a method like int bar(const double&)
  static_assert(has_method<E, int, bar_method_t, const double&>, "failed");
  // E does not have method bar returning double instead of int
  static_assert(!has_method<E, double, bar_method_t, const double&>, "failed");
  // E does not have method taking non-ref const double argument
  static_assert(!has_method<E, int, bar_method_t, const double>, "failed");
  // E does not have method taking non-const double ref argument
  static_assert(!has_method<E, int, bar_method_t, double&>, "failed");
  // E does not have method taking plain double argument
  static_assert(!has_method<E, int, bar_method_t, double>, "failed");
  // E does not have method const method with correct signature otherwise
  // This test ensures a non-const method does not qualify for a
  // check for a method on the const type
  static_assert(!has_method<const E, int, bar_method_t, const double&>,
                "failed");
  // E does not have a foo method
  static_assert(!has_method<E, int, foo_method_t, const double&>, "failed");

  // E2 doesn't have method like int bar()
  static_assert(!has_method<E2, int, bar_method_t>, "failed");
  // E2 does not have non-const method with signature int bar(const double&)
  // This means that a const method won't fulfill a non-const method
  // requirement
  static_assert(!has_method<E2, int, bar_method_t, const double&>, "failed");
  // E2 has method of signature int bar(const double&) const
  static_assert(has_method<const E2, int, bar_method_t, const double&>,
                "failed");
  // E2 does not have method taking non-ref const double
  static_assert(!has_method<E2, int, bar_method_t, const double>, "failed");
  // E2 does not have method taking non-const ref double
  static_assert(!has_method<E2, int, bar_method_t, double&>, "failed");
  // E2 does not have method taking plain double
  static_assert(!has_method<E2, int, bar_method_t, double>, "failed");
  // E2 does not have method with char return type
  static_assert(!has_method<const E2, char, bar_method_t, const double&>,
                "failed");
  // E2 does not have foo method
  static_assert(!has_method<const E2, int, foo_method_t, const double&>,
                "failed");

  // E3 does have a method like int bar(const double&) but is private
  static_assert(!has_method<E3, int, bar_method_t, const double&>, "failed");
}

// trait for member named "member_a"
template <typename T>
using member_a_t = decltype(std::declval<T>().member_a);
// trait for member named "member_b"
template <typename T>
using member_b_t = decltype(std::declval<T>().member_b);

struct M {
  int member_a;
  double member_b;
};

struct M2 {
  double member_a;
};

struct M3 {
  char member_a;
};

struct M4 {
  char member_b;
};

class M5 {
  char member_b;
};

BOOST_AUTO_TEST_CASE(TypeTraitsMember) {
  static_assert(has_member<M, member_a_t, int>, "!");
  static_assert(has_member<M, member_b_t, double>, "!");
  // incorrect type
  static_assert(!has_member<M, member_b_t, int>, "!");
  static_assert(!has_member<M, member_a_t, double>, "!");

  static_assert(has_member<M2, member_a_t, double>, "!");
  static_assert(!has_member<M2, member_a_t, int>, "!");

  static_assert(exists<member_a_t, M>, "!");
  static_assert(exists<member_a_t, M2>, "!");
  static_assert(exists<member_a_t, M3>, "!");
  static_assert(!exists<member_a_t, M4>, "!");

  // private member is not detected
  static_assert(!has_member<M5, member_b_t, char>, "!");

  // private member is not detected.
  static_assert(!exists<member_b_t, M5>, "!");
}

template <typename T>
using nested_a_t = typename T::NestedA;
template <typename T>
using nested_b_t = typename T::NestedB;

struct N {
  struct NestedA;
  class NestedB;
};

struct N2 {
  struct NestedA;
};

struct N3 {
  class NestedB;
};

BOOST_AUTO_TEST_CASE(TypeTraitsNestedType) {
  static_assert(exists<nested_a_t, N>, "!");
  static_assert(exists<nested_b_t, N>, "!");

  static_assert(exists<nested_a_t, N2>, "!");
  static_assert(!exists<nested_b_t, N2>, "!");

  static_assert(!exists<nested_a_t, N3>, "!");
  static_assert(exists<nested_b_t, N3>, "!");
}

// trait for member named "member"
template <typename T>
using member_t = decltype(std::declval<T>().member);

// trait for nested type called "Nested"
template <typename T>
using nested_t = typename T::Nested;

// trait for contained template "meta" with two template params
template <typename T>
using meta_t = typename T::template meta<void, void>;

// combine it into a concept
template <typename T>
constexpr bool SomeConcept =
    require<has_method<T, double, foo_method_t, double, int>,
            has_method<const T, bool, bar_method_t, double&&>,
            has_member<T, member_t, bool>, exists<nested_t, T>,
            exists<meta_t, T>>;

struct A {
  bool member;

  struct Nested {};

  template <typename U, typename V>
  struct meta {};

  double foo(double /*unused*/, int /*unused*/) { return 5; }

  bool bar(double&& /*unused*/) const { return true; }
};

struct A2 {
  bool member;

  struct Nested {};

  template <typename U>
  struct meta {};

  double foo(double /*unused*/, int /*unused*/) { return 5; }

  bool bar(double&& /*unused*/) const { return true; }
};

struct B {
  bool different;

  int foo(double /*unused*/) { return 5; }
};

struct C {
  double foo(int /*unused*/) { return 5; }
};

struct D {
  double bar(double /*unused*/) { return 5; }
};

BOOST_AUTO_TEST_CASE(TypeTraitsConcepts) {
  static_assert(SomeConcept<A>, "A does not fulfill \"SomeConcept\"");
  static_assert(!SomeConcept<A2>, "A2 does not fulfill \"SomeConcept\"");
  static_assert(!SomeConcept<B>, "B does fulfill \"SomeConcept\"");
  static_assert(!SomeConcept<C>, "C does fulfill \"SomeConcept\"");
  static_assert(!SomeConcept<D>, "D does fulfill \"SomeConcept\"");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
