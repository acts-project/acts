// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
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

/// @cond detail
namespace detail {

  /// @brief visitor pattern to extract the surface
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  template <typename parameters_t>
  struct ParametersGetter : public boost::static_visitor<const parameters_t&>
  {
  public:
    explicit ParametersGetter(StateType psType) : sType(psType) {}

    /// @brief call operator for parameters extraction using the boost visitor
    /// pattern
    /// @tparam track_state_t Type of the measurement object (templated)
    /// @param edm The edm object for which the parameters will be extracted
    template <typename track_state_t>
    const parameters_t&
    operator()(const track_state_t& edm) const
    {
      switch (sType) :
        {
        case predicted:
          return edm.predictedState;
        case updated:
          return edm.updatedState;
        default :
          return edm.smoothedState;
        }
    }
  };

  /// @brief get method to be used with the visitor pattern
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename parameters_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  parameters_t
  getParamaters(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm, 
                StateType sType)
  {
    static const ParametersGetter pg = ParametersGetter(sType);
    return boost::apply_visitor(pg, edm);
  }

  /// @brief set method to be used with the visitor pattern
  ///
  /// @tparam parameters_t the parameter type to be retrieved
  ///
  /// @param edm is the boost variant edm object
  template <typename parameters_t, BOOST_VARIANT_ENUM_PARAMS(typename T)>
  void
  setParamaters(const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& edm, 
                parmeters_t pars, 
                StateType sType)
  {
    static const ParametersSetter ps = ParametersSetter(pars, sType);
    return boost::apply_visitor(ps, edm);
  }

}  // namespace detail
}  // namespace Acts