// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/TypeList.hpp"

#include <type_traits>
#include <typeinfo>

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(Utilities)

BOOST_AUTO_TEST_CASE(TypeListCreation) {
  struct A {};

  using MyList = Acts::detail::TypeList<float, int, double, A>;

  bool frontIsFloat = std::is_same_v<Acts::detail::front<MyList>, float>;
  BOOST_CHECK(frontIsFloat);

  bool backIsA = std::is_same_v<Acts::detail::back<MyList>, A>;
  BOOST_CHECK(backIsA);

  BOOST_CHECK_EQUAL(Acts::detail::size<MyList>, 4);
}

BOOST_AUTO_TEST_CASE(TypeListPushFront) {
  class A {};
  class B {};
  class C {};

  using BcList = Acts::detail::TypeList<B, C>;
  auto abc = Acts::detail::push_front<BcList, A>{};
  bool frontIsA = std::is_same_v<Acts::detail::front<decltype(abc)>, A>;
  BOOST_CHECK(frontIsA);
}

BOOST_AUTO_TEST_CASE(TypeListPushBack) {
  class A {};
  class B {};
  class C {};

  using AbList = Acts::detail::TypeList<A, B>;
  auto abc = Acts::detail::push_back<AbList, C>{};
  bool backIsC = std::is_same_v<Acts::detail::back<decltype(abc)>, C>;
  BOOST_CHECK(backIsC);
}

template <typename Head, typename... Tail>
void printTypes(
    [[maybe_unused]] const Acts::detail::TypeList<Head, Tail...>& t) {
  std::cout << typeid(Head).name() << '\n';
  if constexpr (sizeof...(Tail) > 0) {
    Acts::detail::TypeList<Tail...> remainingTypes;
    printTypes(remainingTypes);
  }
}

BOOST_AUTO_TEST_CASE(TypeListPrintType) {
  class A {};
  class B {};
  class C {};

  using AbcList = Acts::detail::TypeList<A, B, C>;
  AbcList a{};
  printTypes(a);
}

BOOST_AUTO_TEST_SUITE_END()
