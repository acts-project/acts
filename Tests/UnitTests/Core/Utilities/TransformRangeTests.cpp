// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/TransformRange.hpp"

using namespace Acts;

BOOST_AUTO_TEST_SUITE(TransformRangeTests)

auto checkSameAddresses = [](const auto& orig, auto& act) {
  BOOST_CHECK_EQUAL(orig.size(), act.size());
  auto it = act.begin();
  for (std::size_t i = 0; i < orig.size(); ++i) {
    BOOST_CHECK_EQUAL(orig.at(i).get(), &act.at(i));
    BOOST_CHECK_EQUAL(orig.at(i).get(), &act[i]);
    BOOST_CHECK_EQUAL(*orig.at(i), act.at(i));

    BOOST_CHECK_EQUAL(orig.at(i).get(), &*it);
    BOOST_CHECK_EQUAL(*orig.at(i).get(), *it);

    ++it;
  }

  BOOST_CHECK(it == act.end());
};

BOOST_AUTO_TEST_CASE(TransformRangeDeref) {
  {
    std::vector<std::unique_ptr<int>> v;
    v.push_back(std::make_unique<int>(1));
    v.push_back(std::make_unique<int>(2));
    v.push_back(std::make_unique<int>(3));

    {
      auto r = detail::make_transform_range<detail::Dereference>(v);
      static_assert(std::is_same_v<decltype(r)::value_type, int&>);
      static_assert(std::is_same_v<decltype(r.at(0)), int&>);
      static_assert(std::is_same_v<decltype(r[0]), int&>);
      static_assert(std::is_same_v<decltype(*r.begin()), int&>);

      r.at(2) = 4;
      BOOST_CHECK_EQUAL(*v.at(2), 4);

      checkSameAddresses(v, r);

      const auto& r_cref = r;
      static_assert(std::is_same_v<decltype(r_cref.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(r_cref[0]), const int&>);
      static_assert(std::is_same_v<decltype(*r_cref.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++r_cref.begin())), const int&>);
      checkSameAddresses(v, r_cref);

      auto cr = detail::make_transform_range<detail::ConstDereference>(v);
      static_assert(std::is_same_v<decltype(cr)::value_type, const int&>);
      static_assert(std::is_same_v<decltype(cr.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(cr[0]), const int&>);
      static_assert(std::is_same_v<decltype(*cr.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++cr.begin())), const int&>);
      checkSameAddresses(v, cr);

      const auto& cr_cref = cr;
      static_assert(std::is_same_v<decltype(cr_cref.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(cr_cref[0]), const int&>);
      static_assert(std::is_same_v<decltype(*cr_cref.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++cr_cref.begin())), const int&>);
      checkSameAddresses(v, cr_cref);
    }

    std::vector<int*> raw_v;
    for (auto& ptr : v) {
      raw_v.push_back(ptr.get());
    }
    {
      auto r = detail::make_transform_range<detail::Dereference>(raw_v);
      static_assert(std::is_same_v<decltype(r)::value_type, int&>);
      static_assert(std::is_same_v<decltype(r.at(0)), int&>);
      static_assert(std::is_same_v<decltype(r[0]), int&>);
      static_assert(std::is_same_v<decltype(*r.begin()), int&>);
      static_assert(std::is_same_v<decltype(*(++r.begin())), int&>);
      checkSameAddresses(v, r);

      auto cr = detail::make_transform_range<detail::ConstDereference>(raw_v);
      static_assert(std::is_same_v<decltype(cr)::value_type, const int&>);
      static_assert(std::is_same_v<decltype(cr.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(cr[0]), const int&>);
      static_assert(std::is_same_v<decltype(*cr.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++cr.begin())), const int&>);
      checkSameAddresses(v, r);
    }
  }

  {
    std::vector<std::unique_ptr<const int>> v;
    v.push_back(std::make_unique<const int>(1));
    v.push_back(std::make_unique<const int>(2));
    v.push_back(std::make_unique<const int>(3));

    {
      auto r = detail::make_transform_range<detail::Dereference>(v);
      static_assert(std::is_same_v<decltype(r)::value_type, const int&>);
      static_assert(std::is_same_v<decltype(r.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(r[0]), const int&>);
      static_assert(std::is_same_v<decltype(*r.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++r.begin())), const int&>);
      checkSameAddresses(v, r);
    }

    std::vector<const int*> raw_v;
    for (auto& ptr : v) {
      raw_v.push_back(ptr.get());
    }
    {
      auto r = detail::make_transform_range<detail::Dereference>(raw_v);
      static_assert(std::is_same_v<decltype(r)::value_type, const int&>);
      static_assert(std::is_same_v<decltype(r.at(0)), const int&>);
      static_assert(std::is_same_v<decltype(r[0]), const int&>);
      static_assert(std::is_same_v<decltype(*r.begin()), const int&>);
      static_assert(std::is_same_v<decltype(*(++r.begin())), const int&>);
      checkSameAddresses(v, r);
    }
  }
}

BOOST_AUTO_TEST_CASE(TransformRangeFromConstRef) {
  std::vector<std::unique_ptr<int>> v;
  v.push_back(std::make_unique<int>(1));
  v.push_back(std::make_unique<int>(2));
  v.push_back(std::make_unique<int>(3));

  const auto& cv = v;

  {
    auto r_cv = detail::make_transform_range<detail::Dereference>(cv);
    static_assert(std::is_same_v<decltype(r_cv)::value_type, const int&>);
    static_assert(std::is_same_v<decltype(r_cv.at(0)), const int&>);
    static_assert(std::is_same_v<decltype(r_cv[0]), const int&>);
    static_assert(std::is_same_v<decltype(*r_cv.begin()), const int&>);
    static_assert(std::is_same_v<decltype(*(++r_cv.begin())), const int&>);

    checkSameAddresses(v, r_cv);
    checkSameAddresses(cv, r_cv);
  }

  {
    auto r_cv = detail::make_transform_range<detail::ConstDereference>(cv);
    static_assert(std::is_same_v<decltype(r_cv)::value_type, const int&>);
    static_assert(std::is_same_v<decltype(r_cv.at(0)), const int&>);
    static_assert(std::is_same_v<decltype(r_cv[0]), const int&>);
    static_assert(std::is_same_v<decltype(*r_cv.begin()), const int&>);
    static_assert(std::is_same_v<decltype(*(++r_cv.begin())), const int&>);

    checkSameAddresses(v, r_cv);
    checkSameAddresses(cv, r_cv);
  }
}

BOOST_AUTO_TEST_SUITE_END()
