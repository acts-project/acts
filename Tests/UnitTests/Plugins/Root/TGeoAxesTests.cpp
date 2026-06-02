// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsPlugins/Root/TGeoAxes.hpp"

#include <stdexcept>
#include <string_view>

using ActsPlugins::TGeoAxes;

// Compile-time construction of valid axes — if these declarations compile,
// the consteval validation accepted them.
constexpr TGeoAxes kXYZ{"XYZ"};
constexpr TGeoAxes kXZY{"XZY"};
constexpr TGeoAxes kYXZ{"YXZ"};
constexpr TGeoAxes kYZX{"YZX"};
constexpr TGeoAxes kZXY{"ZXY"};
constexpr TGeoAxes kZYX{"ZYX"};
// Lowercase (negative direction) variants
constexpr TGeoAxes kxYZ{"xYZ"};
constexpr TGeoAxes kXyZ{"XyZ"};
constexpr TGeoAxes kXYz{"XYz"};
constexpr TGeoAxes kxyz{"xyz"};
constexpr TGeoAxes kXzY{"XzY"};

// The following lines must NOT compile. Uncomment to verify:
// repeated axis
// constexpr TGeoAxes bad_repeated{"XXY"};
// same axis, different case
// constexpr TGeoAxes bad_repeated2{"xXY"};
// invalid character
// constexpr TGeoAxes bad_char{"ABZ"};
// digit not allowed
// constexpr TGeoAxes bad_char2{"X1Z"};
// type error: const char[3] != const char[4]
// constexpr TGeoAxes bad_short{"XY"};
// type error: const char[5] != const char[4]
// constexpr TGeoAxes bad_long{"XYZW"};

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(RootSuite)

BOOST_AUTO_TEST_CASE(CompileTime_ValidLiterals) {
  // Verify value() round-trips correctly for the module-level constexpr
  // objects.
  BOOST_CHECK_EQUAL(kXYZ.value(), "XYZ");
  BOOST_CHECK_EQUAL(kXZY.value(), "XZY");
  BOOST_CHECK_EQUAL(kYXZ.value(), "YXZ");
  BOOST_CHECK_EQUAL(kYZX.value(), "YZX");
  BOOST_CHECK_EQUAL(kZXY.value(), "ZXY");
  BOOST_CHECK_EQUAL(kZYX.value(), "ZYX");
}

BOOST_AUTO_TEST_CASE(CompileTime_CaseVariants) {
  // Lowercase letters (negative direction) are valid.
  BOOST_CHECK_EQUAL(kxYZ.value(), "xYZ");
  BOOST_CHECK_EQUAL(kXyZ.value(), "XyZ");
  BOOST_CHECK_EQUAL(kXYz.value(), "XYz");
  BOOST_CHECK_EQUAL(kxyz.value(), "xyz");
  BOOST_CHECK_EQUAL(kXzY.value(), "XzY");
}

BOOST_AUTO_TEST_CASE(Parse_Valid) {
  // All 6 axis orderings should parse successfully.
  BOOST_CHECK_EQUAL(TGeoAxes::parse("XYZ").value(), "XYZ");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("XZY").value(), "XZY");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("YXZ").value(), "YXZ");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("YZX").value(), "YZX");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("ZXY").value(), "ZXY");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("ZYX").value(), "ZYX");
  // Mixed case
  BOOST_CHECK_EQUAL(TGeoAxes::parse("xYZ").value(), "xYZ");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("XzY").value(), "XzY");
  BOOST_CHECK_EQUAL(TGeoAxes::parse("xyz").value(), "xyz");
}

BOOST_AUTO_TEST_CASE(Parse_WrongLength) {
  BOOST_CHECK_THROW(TGeoAxes::parse("XY"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("XYZW"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse(""), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("X"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Parse_BadChars) {
  BOOST_CHECK_THROW(TGeoAxes::parse("ABZ"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("X1Z"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("XY "),
                    std::invalid_argument);  // trailing space
  BOOST_CHECK_THROW(TGeoAxes::parse("XYW"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Parse_RepeatedAxis) {
  // Same axis letter repeated
  BOOST_CHECK_THROW(TGeoAxes::parse("XXY"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("XYY"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("ZZX"), std::invalid_argument);
  // Same axis, different case (x and X are the same axis)
  BOOST_CHECK_THROW(TGeoAxes::parse("xXY"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("XyY"), std::invalid_argument);
  BOOST_CHECK_THROW(TGeoAxes::parse("zZX"), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
