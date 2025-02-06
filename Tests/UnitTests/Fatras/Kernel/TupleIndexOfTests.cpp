// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Kernel/InteractionList.hpp"

#include <tuple>

namespace ActsFatras::detail {
template <class T, class Tuple>
struct TupleIndexOf;
}  // namespace ActsFatras::detail

using ActsFatras::detail::TupleIndexOf;

BOOST_AUTO_TEST_SUITE(FatrasTupleIndexOf)

BOOST_AUTO_TEST_CASE(Regular) {
  using T = std::tuple<int, double, float>;

  BOOST_CHECK_EQUAL((TupleIndexOf<int, T>::value), 0u);
  BOOST_CHECK_EQUAL((TupleIndexOf<double, T>::value), 1u);
  BOOST_CHECK_EQUAL((TupleIndexOf<float, T>::value), 2u);
}

BOOST_AUTO_TEST_CASE(Duplicates) {
  using T = std::tuple<double, int, double>;

  // should return the first matching type
  BOOST_CHECK_EQUAL((TupleIndexOf<double, T>::value), 0u);
  BOOST_CHECK_EQUAL((TupleIndexOf<int, T>::value), 1u);
}

BOOST_AUTO_TEST_SUITE_END()
