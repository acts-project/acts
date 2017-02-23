// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_DEFAULTPARAMETERDEFINITIONS_H
#define ACTS_DEFAULTPARAMETERDEFINITIONS_H 1

// STL include(s)
#include <cmath>

// ACTS includes
#include "ACTS/Utilities/ParameterTypes.hpp"

namespace Acts {
enum ParDef : unsigned int {
  eLOC_0    = 0,  ///< first coordinate in local surface frame
  eLOC_1    = 1,  ///< second coordinate in local surface frame
  eLOC_R    = eLOC_0,
  eLOC_PHI  = eLOC_1,
  eLOC_RPHI = eLOC_0,
  eLOC_Z    = eLOC_1,
  eLOC_X    = eLOC_0,
  eLOC_Y    = eLOC_1,
  eLOC_D0   = eLOC_0,
  eLOC_Z0   = eLOC_1,
  ePHI      = 2,  ///< phi direction of momentum in global frame
  eTHETA    = 3,  ///< theta direction of momentum in global frame
  eQOP = 4,  ///< charge/momentum for charged tracks, for neutral tracks it is
             /// 1/momentum
  NGlobalPars
};

typedef ParDef ParID_t;
typedef double ParValue_t;

template <ParID_t>
struct par_type;

template <ParID_t par>
using par_type_t = typename par_type<par>::type;

template <>
struct par_type<ParDef::eLOC_0>
{
  typedef local_parameter type;
};

template <>
struct par_type<ParDef::eLOC_1>
{
  typedef local_parameter type;
};

template <>
struct par_type<ParDef::ePHI>
{
  static constexpr double
  pMin()
  {
    return -M_PI;
  }
  static constexpr double
  pMax()
  {
    return M_PI;
  }
  typedef cyclic_parameter<double, pMin, pMax> type;
};

template <>
struct par_type<ParDef::eTHETA>
{
  static constexpr double
  pMin()
  {
    return 0;
  }
  static constexpr double
  pMax()
  {
    return M_PI;
  }
  typedef bound_parameter<double, pMin, pMax> type;
};

template <>
struct par_type<ParDef::eQOP>
{
  typedef unbound_parameter type;
};
}  // end of namespace Acts

#endif  // ACTS_DEFAULTPARAMETERDEFINITIONS_H
