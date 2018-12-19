// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// boost include(s)
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
// forward declaration
template <typename identifier_t,
          typename parameters_t,
          typename jacobian_t,
          ParID_t... params>
class TrackState;

/// @cond detail
namespace detail {
  ///
  /// @brief generate boost::variant type for all possible measurements
  ///
  template <typename identifier_t, typename parameters_t, typename jacobian_t>
  struct trackstate_type_generator;

  /// @cond
  template <typename identifier_t, typename parameters_t, typename jacobian_t>
  struct trackstate_type_generator
  {
    template <typename... T>
    struct container
    {
    };

    template <typename T, typename U>
    struct add_prepended;

    template <ParID_t first, typename... others>
    struct add_prepended<TrackState<identifier_t,
                                    parameters_t,
                                    jacobian_t,
                                    first>,
                         container<others...>>
    {
      using type = container<typename add_prepended<TrackState<identifier_t,
                                                               parameters_t,
                                                               jacobian_t,
                                                               first>,
                                                    others>::type...,
                             others...>;
    };

    template <ParID_t first, ParID_t... others>
    struct add_prepended<TrackState<identifier_t,
                                    parameters_t,
                                    jacobian_t,
                                    first>,
                         TrackState<identifier_t,
                                    parameters_t,
                                    jacobian_t,
                                    others...>>
    {
      using type = TrackState<identifier_t,
                              parameters_t,
                              jacobian_t,
                              first,
                              others...>;
    };

    template <ParID_t... first>
    struct add_prepended<TrackState<identifier_t,
                                    parameters_t,
                                    jacobian_t,
                                    first...>,
                         boost::mpl::na>
    {
      using type = TrackState<identifier_t, parameters_t, jacobian_t, first...>;
    };

    template <typename T, typename C>
    struct add_to_container;

    template <typename T, typename... others>
    struct add_to_container<T, container<others...>>
    {
      using type = container<T, others...>;
    };

    template <typename T>
    struct generator_impl;

    template <ParID_t first, ParID_t... others>
    struct generator_impl<container<TrackState<identifier_t,
                                               parameters_t,
                                               jacobian_t,
                                               first>,
                                    TrackState<identifier_t,
                                               parameters_t,
                                               jacobian_t,
                                               others...>>>
    {
      using type
          = container<TrackState<identifier_t, parameters_t, jacobian_t, first>,
                      TrackState<identifier_t,
                                 parameters_t,
                                 jacobian_t,
                                 others...>,
                      TrackState<identifier_t,
                                 parameters_t,
                                 jacobian_t,
                                 first,
                                 others...>>;
    };

    template <ParID_t first, typename next, typename... others>
    struct generator_impl<container<TrackState<identifier_t,
                                               parameters_t,
                                               jacobian_t,
                                               first>,
                                    next,
                                    others...>>
    {
      using others_combined =
          typename generator_impl<container<next, others...>>::type;
      using prepended = typename add_prepended<TrackState<identifier_t,
                                                          parameters_t,
                                                          jacobian_t,
                                                          first>,
                                               others_combined>::type;
      using type = typename add_to_container<TrackState<identifier_t,
                                                        parameters_t,
                                                        jacobian_t,
                                                        first>,
                                             prepended>::type;
    };

    template <ParID_t v, typename C>
    struct add_to_value_container;

    template <ParID_t v, ParID_t... others>
    struct add_to_value_container<v, std::integer_sequence<ParID_t, others...>>
    {
      using type = std::integer_sequence<ParID_t, others..., v>;
    };

    template <typename T, unsigned int N>
    struct tparam_generator
    {
      using type =
          typename add_to_value_container<static_cast<ParID_t>(N),
                                          typename tparam_generator<T, N - 1>::
                                              type>::type;
    };

    template <typename T>
    struct tparam_generator<T, 0>
    {
      using type = std::integer_sequence<T, static_cast<T>(0)>;
    };

    template <typename T>
    struct converter;

    template <ParID_t... values>
    struct converter<std::integer_sequence<ParID_t, values...>>
    {
      using type = container<TrackState<identifier_t,
                                        parameters_t,
                                        jacobian_t,
                                        values>...>;
    };

    template <typename... types>
    struct to_boost_vector;

    template <typename first, typename... rest>
    struct to_boost_vector<first, rest...>
    {
      using type = typename boost::mpl::
          push_front<typename to_boost_vector<rest...>::type, first>::type;
    };

    template <typename last>
    struct to_boost_vector<last>
    {
      using type = boost::mpl::vector<last>;
    };

    template <typename... MeasTypes>
    struct converter<container<MeasTypes...>>
    {
      using type = typename boost::make_variant_over<
          typename to_boost_vector<MeasTypes...>::type>::type;
    };

    using par_list =
        typename tparam_generator<ParID_t, Acts::NGlobalPars - 1>::type;
    using meas_list    = typename converter<par_list>::type;
    using permutations = typename generator_impl<meas_list>::type;
    using type         = typename converter<permutations>::type;
  };
  /// @endcond
}  // namespace details
/// @endcond
}  // namespace Acts
