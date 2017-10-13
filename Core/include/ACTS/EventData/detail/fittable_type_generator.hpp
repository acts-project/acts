// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_FITTABLE_TYPE_GENERATOR_H
#define ACTS_FITTABLE_TYPE_GENERATOR_H 1

// boost include(s)
#include <boost/mpl/vector.hpp>
#include <boost/variant.hpp>

// ACTS include(s)
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {
// forward declaration
template <typename Identifier, ParID_t... params>
class Measurement;

/// @cond detail
namespace detail {
  /**
   * @brief generate boost::variant type for all possible Measurement's
   */
  template <typename ID>
  struct fittable_type_generator;

  /// @cond
  template <typename ID>
  struct fittable_type_generator
  {
    template <typename... T>
    struct container
    {
    };

    template <typename T, typename U>
    struct add_prepended;

    template <ParID_t first, typename... others>
    struct add_prepended<Measurement<ID, first>, container<others...>>
    {
      typedef container<
          typename add_prepended<Measurement<ID, first>, others>::type...,
          others...>
          type;
    };

    template <ParID_t first, ParID_t... others>
    struct add_prepended<Measurement<ID, first>, Measurement<ID, others...>>
    {
      typedef Measurement<ID, first, others...> type;
    };

    template <ParID_t... first>
    struct add_prepended<Measurement<ID, first...>, boost::mpl::na>
    {
      typedef Measurement<ID, first...> type;
    };

    template <typename T, typename C>
    struct add_to_container;

    template <typename T, typename... others>
    struct add_to_container<T, container<others...>>
    {
      typedef container<T, others...> type;
    };

    template <typename T>
    struct generator_impl;

    template <ParID_t first, ParID_t... others>
    struct generator_impl<container<Measurement<ID, first>,
                                    Measurement<ID, others...>>>
    {
      typedef container<Measurement<ID, first>,
                        Measurement<ID, others...>,
                        Measurement<ID, first, others...>>
          type;
    };

    template <ParID_t first, typename next, typename... others>
    struct generator_impl<container<Measurement<ID, first>, next, others...>>
    {
      typedef typename generator_impl<container<next, others...>>::type
          others_combined;
      typedef
          typename add_prepended<Measurement<ID, first>, others_combined>::type
              prepended;
      typedef typename add_to_container<Measurement<ID, first>, prepended>::type
          type;
    };

    template <ParID_t v, typename C>
    struct add_to_value_container;

    template <ParID_t v, ParID_t... others>
    struct add_to_value_container<v, std::integer_sequence<ParID_t, others...>>
    {
      typedef std::integer_sequence<ParID_t, others..., v> type;
    };

    template <typename T, unsigned int N>
    struct tparam_generator
    {
      typedef
          typename add_to_value_container<static_cast<ParID_t>(N),
                                          typename tparam_generator<T, N - 1>::
                                              type>::type type;
    };

    template <typename T>
    struct tparam_generator<T, 0>
    {
      typedef std::integer_sequence<T, static_cast<T>(0)> type;
    };

    template <typename T>
    struct converter;

    template <ParID_t... values>
    struct converter<std::integer_sequence<ParID_t, values...>>
    {
      typedef container<Measurement<ID, values>...> type;
    };

    template <typename... types>
    struct to_boost_vector;

    template <typename first, typename... rest>
    struct to_boost_vector<first, rest...>
    {
      typedef typename boost::mpl::
          push_front<typename to_boost_vector<rest...>::type, first>::type type;
    };

    template <typename last>
    struct to_boost_vector<last>
    {
      typedef boost::mpl::vector<last> type;
    };

    template <typename... MeasTypes>
    struct converter<container<MeasTypes...>>
    {
      typedef typename boost::make_variant_over<
          typename to_boost_vector<MeasTypes...>::type>::type type;
    };

    typedef typename tparam_generator<ParID_t, Acts::NGlobalPars - 1>::type
                                                     par_list;
    typedef typename converter<par_list>::type       meas_list;
    typedef typename generator_impl<meas_list>::type permutations;
    typedef typename converter<permutations>::type   type;
  };
  /// @endcond
}  // end of namespace details
/// @endcond
}  // end of namespace Acts

#endif  // ACTS_FITTABLE_TYPE_GENERATOR_H
