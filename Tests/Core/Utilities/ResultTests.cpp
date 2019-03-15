// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Rectangle Bounds Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include "Acts/Utilities/Result.hpp"

#include <iostream>
#include <string>
#include <system_error>

using namespace std::string_literals;

namespace {
enum class MyError { Failure = 1, SomethingElse = 2 };

std::error_code
make_error_code(MyError e)
{
  return {static_cast<int>(e), std::generic_category()};
}
}

namespace std {
// register with STL
template <>
struct is_error_code_enum<MyError> : std::true_type
{
};
}

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_SUITE(Surfaces)

  BOOST_AUTO_TEST_CASE(TestConstruction)
  {

    {
      using Result = Result<int, char>;

      Result res = Result::success(42);
      BOOST_CHECK_EQUAL(*res, 42);
      BOOST_CHECK_EQUAL(res.unwrap(), 42);
      BOOST_CHECK(res.ok());
      res = Result::success('e');
      BOOST_CHECK_EQUAL(*res, 'e');
      BOOST_CHECK_EQUAL(res.unwrap(), 'e');
      BOOST_CHECK(res.ok());
      res = Result::failure(42);
      BOOST_CHECK(!res.ok());
      BOOST_CHECK_EQUAL(res.error(), 42);
      BOOST_CHECK_THROW(res.unwrap(), std::runtime_error);
      res = Result::failure('e');
      BOOST_CHECK(!res.ok());
      BOOST_CHECK_EQUAL(res.error(), 'e');
      BOOST_CHECK_THROW(res.unwrap(), std::runtime_error);
    }

    {
      using Result = Result<double, std::string>;

      Result res1("hallo");
      BOOST_CHECK(!res1.ok());
      BOOST_CHECK_EQUAL(res1.error(), "hallo");
      BOOST_CHECK_THROW(res1.unwrap(), std::runtime_error);
      res1 = Result::failure("hallo");
      BOOST_CHECK(!res1.ok());
      BOOST_CHECK_EQUAL(res1.error(), "hallo");
      BOOST_CHECK_THROW(res1.unwrap(), std::runtime_error);

      Result res2(4.5);
      BOOST_CHECK(res2.ok());
      BOOST_CHECK(*res2 == 4.5);
      BOOST_CHECK(res2.unwrap() == 4.5);
      res2 = Result::success(4.5);
      BOOST_CHECK(res2.ok());
      BOOST_CHECK_EQUAL(*res2, 4.5);
      BOOST_CHECK_EQUAL(res2.unwrap(), 4.5);
    }
  }

  BOOST_AUTO_TEST_CASE(TestErrorCodes)
  {
    auto err1 = MyError::Failure;

    std::error_code ec = err1;

    {
      using Result = Result<double, MyError>;

      Result res(42);
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, 42.);

      Result res2(err1);
      BOOST_CHECK(!res2.ok());
      BOOST_CHECK_EQUAL(res2.error(), err1);
      BOOST_CHECK_THROW(res2.unwrap(), std::runtime_error);
    }

    {
      using Result = Result<double, std::error_code>;

      Result res(42);
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, 42.);
      BOOST_CHECK_EQUAL(res.unwrap(), 42);
      res = 46;
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, 46.);
      BOOST_CHECK_EQUAL(res.unwrap(), 46);

      Result res2(ec);
      BOOST_CHECK(!res2.ok());
      BOOST_CHECK_EQUAL(res2.error(), ec);
      BOOST_CHECK_EQUAL(res2.error(), err1);

      res2 = MyError::SomethingElse;
      BOOST_CHECK(!res2.ok());
      BOOST_CHECK_EQUAL(res2.error(), MyError::SomethingElse);
      BOOST_CHECK(res2.error() != MyError::Failure);
    }

    {
      using Result = Result<double>;

      Result res(42);
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, 42.);
      BOOST_CHECK_EQUAL(res.unwrap(), 42);
      res = 46;
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, 46.);
      BOOST_CHECK_EQUAL(res.unwrap(), 46);

      Result res2(ec);
      BOOST_CHECK(!res2.ok());
      BOOST_CHECK_EQUAL(res2.error(), ec);
      BOOST_CHECK_EQUAL(res2.error(), err1);

      res2 = MyError::SomethingElse;
      BOOST_CHECK(!res2.ok());
      BOOST_CHECK_EQUAL(res2.error(), MyError::SomethingElse);
      BOOST_CHECK(res2.error() != MyError::Failure);
    }

    {
      using Result = Result<std::string>;

      Result res("hallo");

      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, "hallo");
      BOOST_CHECK_EQUAL(res.unwrap(), "hallo");

      res = "something else";
      BOOST_CHECK(res.ok());
      BOOST_CHECK_EQUAL(*res, "something else");
      BOOST_CHECK_EQUAL(res.unwrap(), "something else");

      res = MyError::SomethingElse;
      BOOST_CHECK(!res.ok());
      BOOST_CHECK_EQUAL(res.error(), MyError::SomethingElse);
      BOOST_CHECK(res.error() != MyError::Failure);
    }
  }

  struct NoCopy
  {
    NoCopy(int i) : num(i){};
    NoCopy(const NoCopy&) = delete;
    NoCopy&
    operator=(const NoCopy&)
        = delete;
    NoCopy(NoCopy&&) = default;

    int num;
  };

  BOOST_AUTO_TEST_CASE(CopyBehavious)
  {

    using Result = Result<NoCopy>;

    NoCopy n(5);
    Result res = std::move(n);
    BOOST_CHECK(res.ok());
  }

  BOOST_AUTO_TEST_SUITE_END()
}
}
