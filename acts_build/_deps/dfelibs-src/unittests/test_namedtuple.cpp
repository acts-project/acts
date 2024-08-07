// SPDX-License-Identifier: MIT

/// \file
/// \brief Unit tests for NamedTuple

#include <boost/test/unit_test.hpp>

#include "dfe/dfe_namedtuple.hpp"
#include "record.hpp"

BOOST_TEST_DONT_PRINT_LOG_VALUE(Record::Tuple)

BOOST_AUTO_TEST_CASE(namedtuple_tuple_helpers) {
  using std::tuple_element_t;
  using std::tuple_size;
  using Tuple = Record::Tuple;

  BOOST_TEST(tuple_size<Tuple>::value == 7u);
  BOOST_TEST((std::is_same<tuple_element_t<0, Tuple>, int16_t>::value));
  BOOST_TEST((std::is_same<tuple_element_t<1, Tuple>, int32_t>::value));
  BOOST_TEST((std::is_same<tuple_element_t<2, Tuple>, int64_t>::value));
  BOOST_TEST((std::is_same<tuple_element_t<3, Tuple>, uint64_t>::value));
  BOOST_TEST((std::is_same<tuple_element_t<4, Tuple>, float>::value));
  BOOST_TEST((std::is_same<tuple_element_t<5, Tuple>, double>::value));
  BOOST_TEST((std::is_same<tuple_element_t<6, Tuple>, bool>::value));
}

BOOST_AUTO_TEST_CASE(namedtuple_names) {
  auto example = make_record(123);
  BOOST_TEST(example.names().size() == 7);
  BOOST_TEST(example.names().at(0) == "x");
  BOOST_TEST(example.names().at(1) == "y");
  BOOST_TEST(example.names().at(2) == "z");
  BOOST_TEST(example.names().at(3) == "a");
  BOOST_TEST(example.names().at(4) == "b");
  BOOST_TEST(example.names().at(5) == "c");
  BOOST_TEST(example.names().at(6) == "d");
}

BOOST_AUTO_TEST_CASE(nametuple_assign_from_tuple) {
  // check default values
  Record r;
  BOOST_TEST(r.x == 0);
  BOOST_TEST(r.y == 0);
  BOOST_TEST(r.z == 0);
  BOOST_TEST(r.a == 0);
  BOOST_TEST(r.b == 0.0f);
  BOOST_TEST(r.c == 0.0);
  BOOST_TEST(not r.d);
  // assign NamedTuple from a regular tuple
  r = std::make_tuple(-1, 1, 2, -3, 1.23f, 6.54, true);
  BOOST_TEST(r.x == -1);
  BOOST_TEST(r.y == 1);
  BOOST_TEST(r.z == 2);
  BOOST_TEST(r.a == -3);
  BOOST_TEST(r.b == 1.23f);
  BOOST_TEST(r.c == 6.54);
  BOOST_TEST(r.d);
}

BOOST_AUTO_TEST_CASE(namedtuple_assign_to_tuple) {
  using std::get;

  // assign regular tuple from namedtuple w/ default values
  Record::Tuple t = Record();
  BOOST_TEST(get<0>(t) == 0);
  BOOST_TEST(get<1>(t) == 0);
  BOOST_TEST(get<2>(t) == 0);
  BOOST_TEST(get<3>(t) == 0);
  BOOST_TEST(get<4>(t) == 0.0f);
  BOOST_TEST(get<5>(t) == 0.0);
  BOOST_TEST(not get<6>(t));
  // assign regular tuple from namedtuple
  t = make_record(42);
  Record cmp = make_record(42);
  BOOST_TEST(get<0>(t) == cmp.x);
  BOOST_TEST(get<1>(t) == cmp.y);
  BOOST_TEST(get<2>(t) == cmp.z);
  BOOST_TEST(get<3>(t) == cmp.a);
  BOOST_TEST(get<4>(t) == cmp.b);
  BOOST_TEST(get<5>(t) == cmp.c);
  BOOST_TEST(get<6>(t) == cmp.d);
}

BOOST_AUTO_TEST_CASE(namedtuple_get) {
  using std::get;

  auto r = make_record(42);
  BOOST_TEST(r.x == get<0>(r));
  BOOST_TEST(r.y == get<1>(r));
  BOOST_TEST(r.z == get<2>(r));
  BOOST_TEST(r.a == get<3>(r));
  BOOST_TEST(r.b == get<4>(r));
  BOOST_TEST(r.c == get<5>(r));
  BOOST_TEST(r.d == get<6>(r));
}

BOOST_AUTO_TEST_CASE(namedtuple_get_assign) {
  using std::get;

  auto r = make_record(42);
  BOOST_TEST(r.x == get<0>(r));
  BOOST_TEST(r.x == get<0>(r));
  auto updated = r.x + 2;
  get<0>(r) = updated;
  BOOST_TEST(r.x == updated);
  BOOST_TEST(get<0>(r) == updated);
}
