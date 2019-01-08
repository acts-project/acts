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
  namespace detail {
    struct get_surface_visitor : public boost::static_visitor<const Surface*>
    {
      template <typename measurement_t>
      const Surface*
      operator()(const measurement_t& meas) const
      {
        return &meas.referenceSurface();
      }
    };
  }

  template <typename M>
  const Surface*
  getSurface(const M& meas)
  {
    return boost::apply_visitor(detail::get_surface_visitor(), meas);
  }
}
}
