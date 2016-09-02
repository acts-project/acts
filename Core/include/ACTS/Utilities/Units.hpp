// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_UNITS_HPP
#define ACTS_UNITS_HPP 1

namespace Acts {

/// @brief Unit and conversion constants
///
/// In order to make sure that always consistent numbers and units are used in
/// calculations, one should make use of the constants defined in this
/// namespace to express the units of numerical values. The following
/// conventions are used:
/// - input values to ACTS must be given as `numerical_value * unit_constant`
/// - output values can be converted into the desired unit using
///   `numerical_value / unit_constant`
///
/// Examples:
/// @code
/// using namespace Acts::units;
/// // specify input variables
/// double momentum = 2.5 * _GeV;
/// double width = 23 * _cm;
/// double velocity = 345 * _m/_s;
/// double density = 1.2 * _kg/(_m*_m*_m);
/// double bfield = 2 * _T;
///
/// // convert output values
/// double x_position = trackPars.position().x() / _mm;
/// double pt         = trackPars.momentum().pT() / _TeV;
/// @endcode
namespace units {

/// @name length units
/// @{
#ifdef DOXYGEN
  const double _m = unspecified;
#else
  const double _m   = 1e3;
#endif  // DOXYGEN
  const double _km = 1e3 * _m;
  const double _cm = 1e-2 * _m;
  const double _mm = 1e-3 * _m;
  const double _um = 1e-6 * _m;
  const double _nm = 1e-9 * _m;
  const double _pm = 1e-12 * _m;
  const double _fm = 1e-15 * _m;
/// @}

/// @name mass units
/// @{
#ifdef DOXYGEN
  const double _kg = unspecified;
#else
  const double _kg  = 1e3;
#endif  // DOXYGEN
  const double _g  = 1e-3 * _kg;
  const double _mg = 1e-6 * _kg;
  /// atomic mass unit
  const double _u = 1.660539040e-27 * _kg;
/// @}

/// @name time units
/// @{
#ifdef DOXYGEN
  const double _s = unspecified;
#else
  const double _s   = 1;
#endif  // DOXYGEN
  const double _ms = 1e-3 * _s;
  const double _h  = 3600 * _s;
/// @}

/// @name energy units
/// @{
#ifdef DOXYGEN
  const double _MeV = unspecified;
#else
  const double _MeV = 1e-3;
#endif  // DOXYGEN
  const double _GeV = 1e3 * _MeV;
  const double _TeV = 1e6 * _MeV;
  const double _keV = 1e-3 * _MeV;
  const double _eV  = 1e-6 * _MeV;
/// @}

/// @name charge units
/// @{
#ifdef DOXYGEN
  const double _C = unspecified;
#else
  const double _C   = 1. / 1.60217733e-19;
#endif  // DOXYGEN
  /// elementary electric charge
  const double _e = 1.60217733e-19 * _C;
  /// @}

  /// @name derived units
  /// @{
  const double _N = _kg * _m / (_s * _s);
  const double _J = _N * _m;
  const double _T = _kg / (_C * _s);
  /// @}

  /// @name fundamental physical constants in SI units
  /// @{
  const double _c         = 2.99792458e8 * _m / _s;
  const double _hbar      = 1.05457266e-34 * _J * _s;
  const double _el_charge = _e / _C;
  /// @}

  /// @cond
  /// @brief internally used conversion constants
  namespace {
    // 1 GeV = 1e9 * e * 1 V = 1.60217733e-10 As * 1 V = 1.60217733e-10 J
    // ==> 1 J = 1 / 1.60217733e-10 GeV
    const double _GeV_per_J = _GeV / (_el_charge * 1e9 * _J);
    // hbar * c = 3.161529298809983e-26 J * m
    // ==> hbar * c * _GeV_per_J = 1.973270523563071e-16 GeV * m
    const double _mm_times_GeV = _c * _hbar * _GeV_per_J;
  }
  /// @endcond

  /// @brief physical quantities
  ///
  /// These constants can be used to select the correct conversion function
  /// from SI units to natural units and vice versa.
  enum Quantity { MOMENTUM, ENERGY, LENGTH, MASS };

  template <Quantity>
  double
  SI2Nat(const double);

  template <Quantity>
  double
  Nat2SI(const double);

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
#endif  // ACTS_UNITS_HPP
