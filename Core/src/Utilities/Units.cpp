// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace units {
  template <>
  double
  SI2Nat<ENERGY>(const double E)
  {
    static const double conversion = _GeV_per_J;
    return E * conversion;
  }

  template <>
  double
  Nat2SI<ENERGY>(const double E)
  {
    static const double conversion = 1. / _GeV_per_J;
    return E * conversion;
  }

  template <>
  double
  SI2Nat<LENGTH>(const double l)
  {
    static const double conversion = 1. / _mm_times_GeV;
    return l * conversion;
  }

  template <>
  double
  Nat2SI<LENGTH>(const double l)
  {
    static const double conversion = _mm_times_GeV;
    return l * conversion;
  }

  template <>
  double
  SI2Nat<MOMENTUM>(const double p)
  {
    // p(NU) = p(SI) * c
    static const double conversion = _c * _GeV_per_J;
    return p * conversion;
  }

  template <>
  double
  Nat2SI<MOMENTUM>(const double p)
  {
    // p(SI) = p(NU)/c
    static const double conversion = 1. / (_c * _GeV_per_J);
    return p * conversion;
  }

  template <>
  double
  SI2Nat<MASS>(const double m)
  {
    // p(NU) = p(SI) * c
    static const double conversion = _c * _c * _GeV_per_J;
    return m * conversion;
  }

  template <>
  double
  Nat2SI<MASS>(const double m)
  {
    // p(SI) = p(NU)/c
    static const double conversion = 1. / (_c * _c * _GeV_per_J);
    return m * conversion;
  }
}  // namespace units
}  // namespace Acts
