// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <memory>
#include <random>
#include <variant>

#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/EventData/detail/fittable_type_generator.hpp"

#include <boost/hana/equal.hpp>
#include <boost/hana/integral_constant.hpp>
#include <boost/hana/tuple.hpp>
#include <boost/hana/type.hpp>

namespace hana = boost::hana;
using namespace hana::literals;

#define _T(...) hana::tuple_c<size_t, __VA_ARGS__>

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_CASE(index_combination_generation_test) {
  {
    constexpr auto result = detail::unique_ordered_sublists<1>();
    constexpr auto expected = hana::make_tuple(_T(0));
    static_assert(result == expected, "At size 1 is not equal");
  }  // namespace Test
  {
    constexpr auto result = detail::unique_ordered_sublists<2>();
    constexpr auto expected = hana::make_tuple(_T(0), _T(1), _T(0, 1));
    static_assert(result == expected, "At size 2 is not equal");
  }
  {
    constexpr auto result = detail::unique_ordered_sublists<3>();
    constexpr auto expected = hana::make_tuple(_T(0), _T(1), _T(0, 1), _T(2),
                                               _T(0, 2), _T(1, 2), _T(0, 1, 2));
    static_assert(result == expected, "At size 3 is not equal");
  }
  {
    constexpr auto result = detail::unique_ordered_sublists<4>();
    constexpr auto expected =
        hana::make_tuple(_T(0), _T(1), _T(0, 1), _T(2), _T(0, 2), _T(1, 2),
                         _T(0, 1, 2), _T(3), _T(0, 3), _T(1, 3), _T(0, 1, 3),
                         _T(2, 3), _T(0, 2, 3), _T(1, 2, 3), _T(0, 1, 2, 3));
    static_assert(result == expected, "At size 4 is not equal");
  }
  {
    constexpr auto result = detail::unique_ordered_sublists<5>();
    constexpr auto expected = hana::make_tuple(
        _T(0), _T(1), _T(0, 1), _T(2), _T(0, 2), _T(1, 2), _T(0, 1, 2), _T(3),
        _T(0, 3), _T(1, 3), _T(0, 1, 3), _T(2, 3), _T(0, 2, 3), _T(1, 2, 3),
        _T(0, 1, 2, 3), _T(4), _T(0, 4), _T(1, 4), _T(0, 1, 4), _T(2, 4),
        _T(0, 2, 4), _T(1, 2, 4), _T(0, 1, 2, 4), _T(3, 4), _T(0, 3, 4),
        _T(1, 3, 4), _T(0, 1, 3, 4), _T(2, 3, 4), _T(0, 2, 3, 4),
        _T(1, 2, 3, 4), _T(0, 1, 2, 3, 4));
    static_assert(result == expected, "At size 5 is not equal");
  }
}  // namespace Acts

using par_t = ParID_t;
using Source = MinimalSourceLink;

template <par_t... pars>
struct meas_factory {
  using type = Measurement<Source, par_t, pars...>;
};

constexpr par_t operator"" _p(unsigned long long i) {
  return par_t(i);
}

BOOST_AUTO_TEST_CASE(variant_bound_measurement_generation_test) {
  {
    using actual = detail::type_generator_t<par_t, meas_factory>;
    using expected = std::variant<
        Measurement<Source, par_t, 0_p>, Measurement<Source, par_t, 1_p>,
        Measurement<Source, par_t, 0_p, 1_p>, Measurement<Source, par_t, 2_p>,
        Measurement<Source, par_t, 0_p, 2_p>,
        Measurement<Source, par_t, 1_p, 2_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p>,
        Measurement<Source, par_t, 3_p>, Measurement<Source, par_t, 0_p, 3_p>,
        Measurement<Source, par_t, 1_p, 3_p>,
        Measurement<Source, par_t, 0_p, 1_p, 3_p>,
        Measurement<Source, par_t, 2_p, 3_p>,
        Measurement<Source, par_t, 0_p, 2_p, 3_p>,
        Measurement<Source, par_t, 1_p, 2_p, 3_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 3_p>,
        Measurement<Source, par_t, 4_p>, Measurement<Source, par_t, 0_p, 4_p>,
        Measurement<Source, par_t, 1_p, 4_p>,
        Measurement<Source, par_t, 0_p, 1_p, 4_p>,
        Measurement<Source, par_t, 2_p, 4_p>,
        Measurement<Source, par_t, 0_p, 2_p, 4_p>,
        Measurement<Source, par_t, 1_p, 2_p, 4_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 4_p>,
        Measurement<Source, par_t, 3_p, 4_p>,
        Measurement<Source, par_t, 0_p, 3_p, 4_p>,
        Measurement<Source, par_t, 1_p, 3_p, 4_p>,
        Measurement<Source, par_t, 0_p, 1_p, 3_p, 4_p>,
        Measurement<Source, par_t, 2_p, 3_p, 4_p>,
        Measurement<Source, par_t, 0_p, 2_p, 3_p, 4_p>,
        Measurement<Source, par_t, 1_p, 2_p, 3_p, 4_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 3_p, 4_p>,
        Measurement<Source, par_t, 5_p>, Measurement<Source, par_t, 0_p, 5_p>,
        Measurement<Source, par_t, 1_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 5_p>,
        Measurement<Source, par_t, 2_p, 5_p>,
        Measurement<Source, par_t, 0_p, 2_p, 5_p>,
        Measurement<Source, par_t, 1_p, 2_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 5_p>,
        Measurement<Source, par_t, 3_p, 5_p>,
        Measurement<Source, par_t, 0_p, 3_p, 5_p>,
        Measurement<Source, par_t, 1_p, 3_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 3_p, 5_p>,
        Measurement<Source, par_t, 2_p, 3_p, 5_p>,
        Measurement<Source, par_t, 0_p, 2_p, 3_p, 5_p>,
        Measurement<Source, par_t, 1_p, 2_p, 3_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 3_p, 5_p>,
        Measurement<Source, par_t, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 4_p, 5_p>,
        Measurement<Source, par_t, 1_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 4_p, 5_p>,
        Measurement<Source, par_t, 2_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 2_p, 4_p, 5_p>,
        Measurement<Source, par_t, 1_p, 2_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 4_p, 5_p>,
        Measurement<Source, par_t, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 1_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 2_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 2_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 1_p, 2_p, 3_p, 4_p, 5_p>,
        Measurement<Source, par_t, 0_p, 1_p, 2_p, 3_p, 4_p, 5_p>>;
    static_assert(std::is_same<actual, expected>::value,
                  "Variant is not identical");
  }
}

using freePar_t = FreeParametersIndices;

template <freePar_t... pars>
struct meas_factory2 {
  using type = Measurement<Source, freePar_t, pars...>;
};

constexpr freePar_t operator"" _fp(unsigned long long i) {
  return freePar_t(i);
}

BOOST_AUTO_TEST_CASE(variant_free_measurement_generation_test) {
  {
    using actual = detail::type_generator_t<freePar_t, meas_factory2>;
    using expected = std::variant<
        Measurement<Source, freePar_t, 0_fp>, Measurement<Source, freePar_t, 1_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp>,
        Measurement<Source, freePar_t, 2_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp>,
        Measurement<Source, freePar_t, 3_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp>,
        Measurement<Source, freePar_t, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp>,
        Measurement<Source, freePar_t, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 5_fp>,
        Measurement<Source, freePar_t, 2_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 5_fp>,
        Measurement<Source, freePar_t, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 5_fp>,
        Measurement<Source, freePar_t, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp>,
        Measurement<Source, freePar_t, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 6_fp>,
        Measurement<Source, freePar_t, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 6_fp>,
        Measurement<Source, freePar_t, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 6_fp>,
        Measurement<Source, freePar_t, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp>,
        Measurement<Source, freePar_t, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 7_fp>,
        Measurement<Source, freePar_t, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 7_fp>,
        Measurement<Source, freePar_t, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 7_fp>,
        Measurement<Source, freePar_t, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>,
        Measurement<Source, freePar_t, 0_fp, 1_fp, 2_fp, 3_fp, 4_fp, 5_fp, 6_fp, 7_fp>>;
    static_assert(std::is_same<actual, expected>::value,
                  "Variant is not identical");
  }
}
}  // namespace Test
}  // namespace Acts
