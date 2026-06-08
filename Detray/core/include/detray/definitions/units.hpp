// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"

// System include(s)
#include <cmath>
#include <numbers>

namespace detray {

/// Unit conversion factors
template <concepts::scalar scalar_t>
struct unit {
  /// Length, native unit mm
  /// @{
  static constexpr scalar_t um{static_cast<scalar_t>(1e-3)};
  static constexpr scalar_t mm{static_cast<scalar_t>(1.0)};
  static constexpr scalar_t cm{static_cast<scalar_t>(1e1)};
  static constexpr scalar_t m{static_cast<scalar_t>(1e3)};
  /// @}

  /// Area, native unit mm2
  /// @{
  /// mm²
  static constexpr scalar_t mm2{static_cast<scalar_t>(1.0)};
  /// cm²
  static constexpr scalar_t cm2{static_cast<scalar_t>(1e2)};
  /// m²
  static constexpr scalar_t m2{static_cast<scalar_t>(1e6)};
  /// @}

  /// Volume, native unit mm3
  /// @{
  /// mm³
  static constexpr scalar_t mm3{static_cast<scalar_t>(1.0)};
  /// cm³
  static constexpr scalar_t cm3{static_cast<scalar_t>(1e3)};
  /// m³
  static constexpr scalar_t m3{static_cast<scalar_t>(1e9)};
  /// @}

  /// Time, native unit mm{[speed-of-light * time]{mm/s * s}}
  /// @{
  static constexpr scalar_t s{static_cast<scalar_t>(299792458000.0)};
  static constexpr scalar_t fs{static_cast<scalar_t>(1e-15 * 299792458000.0)};
  static constexpr scalar_t ps{static_cast<scalar_t>(1e-12 * 299792458000.0)};
  static constexpr scalar_t ns{static_cast<scalar_t>(1e-9 * 299792458000.0)};
  static constexpr scalar_t us{static_cast<scalar_t>(1e-6 * 299792458000.0)};
  static constexpr scalar_t ms{static_cast<scalar_t>(1e-3 * 299792458000.0)};
  static constexpr scalar_t min{static_cast<scalar_t>(60.0 * 299792458000.0)};
  static constexpr scalar_t h{static_cast<scalar_t>(3600.0 * 299792458000.0)};
  /// @}

  /// Energy, native unit GeV
  /// @{
  static constexpr scalar_t eV{static_cast<scalar_t>(1e-9)};
  static constexpr scalar_t keV{static_cast<scalar_t>(1e-6)};
  static constexpr scalar_t MeV{static_cast<scalar_t>(1e-3)};
  static constexpr scalar_t GeV{static_cast<scalar_t>(1.0)};
  static constexpr scalar_t TeV{static_cast<scalar_t>(1e3)};
  /// @}

  /// Atomic mass unit u
  /// 1u == 0.93149410242 GeV/c
  static constexpr scalar_t u{static_cast<scalar_t>(0.93149410242)};

  /// Mass
  ///     1eV/c² == 1.782662e-36kg
  ///    1GeV/c² == 1.782662e-27kg
  /// ->     1kg == (1/1.782662e-27)GeV/c²
  /// ->      1g == (1/(1e3*1.782662e-27))GeV/c²
  /// @{
  static constexpr scalar_t g{static_cast<scalar_t>(1.0 / 1.782662e-24)};
  static constexpr scalar_t kg{static_cast<scalar_t>(1.0 / 1.782662e-27)};
  /// @}

  /// Amount of substance, native unit mol
  static constexpr scalar_t mol{static_cast<scalar_t>(1.0)};

  /// Charge, native unit e (elementary charge)
  static constexpr scalar_t e{static_cast<scalar_t>(1.0)};

  /// Magnetic field, native unit GeV/(e*mm)
  static constexpr scalar_t T{static_cast<scalar_t>(
      0.000299792458)};  // equivalent to c in appropriate SI units

  // Angles, native unit radian
  static constexpr scalar_t mrad{static_cast<scalar_t>(1e-3)};
  static constexpr scalar_t rad{static_cast<scalar_t>(1.0)};
  static constexpr scalar_t degree{
      static_cast<scalar_t>(0.017453292519943295)};  // pi / 180
};

/// Physical and mathematical constants
template <concepts::scalar scalar_t>
struct constant {
  /// Euler's number
  static constexpr scalar_t e{static_cast<scalar_t>(std::numbers::e)};
  /// Base 2 logarithm of e
  static constexpr scalar_t log2e{static_cast<scalar_t>(std::numbers::log2e)};
  /// Base 10 logarithm of e
  static constexpr scalar_t log10e{static_cast<scalar_t>(std::numbers::log10e)};
  /// Natural logarithm of 2
  static constexpr scalar_t ln2{static_cast<scalar_t>(std::numbers::ln2)};
  /// Natural logarithm of 10
  static constexpr scalar_t ln10{static_cast<scalar_t>(std::numbers::ln10)};

  /// π
  static constexpr scalar_t pi{static_cast<scalar_t>(std::numbers::pi)};
  /// π/2
  static constexpr scalar_t pi_2{static_cast<scalar_t>(std::numbers::pi / 2.)};
  /// π/4
  static constexpr scalar_t pi_4{static_cast<scalar_t>(std::numbers::pi / 4.)};
  /// 1/π
  static constexpr scalar_t inv_pi{
      static_cast<scalar_t>(1. / std::numbers::pi)};
  /// 2/π
  static constexpr scalar_t inv_2_pi{
      static_cast<scalar_t>(2. / std::numbers::pi)};

  /// √2
  static constexpr scalar_t sqrt2{static_cast<scalar_t>(std::numbers::sqrt2)};
  /// 1/(√2)
  static constexpr scalar_t inv_sqrt2{
      static_cast<scalar_t>(1. / std::numbers::sqrt2)};

  /// Avogadro constant
  static constexpr scalar_t avogadro{
      static_cast<scalar_t>(6.02214076e23 / unit<scalar_t>::mol)};

  /// Reduced Planck constant h/2*pi.
  ///
  /// Computed from CODATA 2018 constants to double precision.
  static constexpr scalar_t hbar{static_cast<scalar_t>(
      6.582119569509066e-25 * unit<scalar_t>::GeV * unit<scalar_t>::s)};

  // electron mass
  // https://physics.nist.gov/cgi-bin/cuu/Value?eqmec2mev (Visited Aug 1st,
  // 2024)
  static constexpr scalar_t m_e{
      static_cast<scalar_t>(0.51099895069 * unit<scalar_t>::MeV)};
};

}  // namespace detray
