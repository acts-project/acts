// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include "Acts/EventData/Measurement.hpp"

namespace Acts {

class Surface;

namespace MeasurementHelpers {

  /// @brief Extract surface from a type erased measurement object
  /// @tparam M The FittableMeasurement type
  /// @return const pointer to the extracted surface
  template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
  const Surface*
  getSurface(
      const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& fittable_measurement)
  {
    return boost::apply_visitor(
        [](const auto& meas) { return &meas.referenceSurface(); },
        fittable_measurement);
  }

  template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
  size_t
  getSize(
      const boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>& fittable_measurement)
  {
    return boost::apply_visitor([](const auto& meas) { return meas.size(); },
                                fittable_measurement);
  }
}
}
