// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <tuple>

namespace Acts::detail {

template <typename cmp_t>
concept ComponentWithoutCovarianceConcept =
    requires { typename std::tuple_size<std::remove_cvref_t<cmp_t>>::type; } &&
    (std::tuple_size_v<std::remove_cvref_t<cmp_t>> == 2) &&
    requires(const cmp_t& cmp) {
      { std::get<0>(cmp) } -> std::convertible_to<const double&>;
      { std::get<1>(cmp) } -> std::convertible_to<const BoundVector&>;
    };

template <typename cmp_t>
concept ComponentWithCovarianceConcept =
    requires { typename std::tuple_size<std::remove_cvref_t<cmp_t>>::type; } &&
    (std::tuple_size_v<std::remove_cvref_t<cmp_t>> == 3) &&
    requires(const cmp_t& cmp) {
      { std::get<0>(cmp) } -> std::convertible_to<const double&>;
      { std::get<1>(cmp) } -> std::convertible_to<const BoundVector&>;
      { std::get<2>(cmp) } -> std::convertible_to<const BoundMatrix&>;
    };

template <typename cmp_t>
concept ComponentConcept = ComponentWithoutCovarianceConcept<cmp_t> ||
                           ComponentWithCovarianceConcept<cmp_t>;

template <typename projector_t, typename cmp_t>
concept ComponentWithoutCovarianceProjectorConcept =
    requires(const projector_t& proj, const cmp_t& cmp) {
      { proj(cmp) } -> ComponentWithoutCovarianceConcept;
    };

template <typename projector_t, typename cmp_t>
concept ComponentWithCovarianceProjectorConcept =
    requires(const projector_t& proj, const cmp_t& cmp) {
      { proj(cmp) } -> ComponentWithCovarianceConcept;
    };

template <typename projector_t, typename cmp_t>
concept ComponentProjectorConcept =
    ComponentWithoutCovarianceProjectorConcept<projector_t, cmp_t> ||
    ComponentWithCovarianceProjectorConcept<projector_t, cmp_t>;

template <typename component_range_t, typename projector_t>
concept ComponentRangeAndProjectorWithoutCovarianceConcept =
    std::ranges::range<component_range_t> &&
    requires(const component_range_t& cmps, const projector_t& proj) {
      { proj(*cmps.begin()) } -> ComponentWithoutCovarianceConcept;
    };

template <typename component_range_t, typename projector_t>
concept ComponentRangeAndProjectorWithCovarianceConcept =
    std::ranges::range<component_range_t> &&
    requires(const component_range_t& cmps, const projector_t& proj) {
      { proj(*cmps.begin()) } -> ComponentWithCovarianceConcept;
    };

}  // namespace Acts::detail
