// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Extendable Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Utilities/detail/Extendable.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // This tests the implementation of the ActionList
  // and the standard aborters
  BOOST_AUTO_TEST_CASE(Extendable_)
  {
    struct type_a
    {
      double va_a = 0.;
    };
    struct type_b
    {
      double va_b = 0.;
    };

    detail::Extendable<> nullist;
    (void)nullist;
    BOOST_TEST(std::tuple_size<std::tuple<>>::value == 0);

    detail::Extendable<type_a> alist;
    auto&                      a0_object = alist.get<type_a>();
    a0_object.va_a                       = 1.;
    BOOST_TEST(alist.get<type_a>().va_a == 1.);

    detail::Extendable<type_a, type_b> ablist;
    auto& a1_object = ablist.get<type_a>();
    a1_object.va_a  = 2.;
    auto& b1_object = ablist.get<type_b>();
    b1_object.va_b  = 3.;
    BOOST_TEST(ablist.get<type_a>().va_a == 2.);
    BOOST_TEST(ablist.get<type_b>().va_b == 3.);
  }

}  // namespace Test
}  // namespace Acts
