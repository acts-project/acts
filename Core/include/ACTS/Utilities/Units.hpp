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
  const double _m   = 1000;
#endif  // DOXYGEN
  const double _km = 1000 * _m;
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
  const double _kg  = 1000;
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
  const double _MeV = 1;
#endif  // DOXYGEN
  const double _GeV = 1e3 * _MeV;
  const double _TeV = 1e6 * _MeV;
  const double _keV = 1e-3 * _MeV;
  const double _eV  = 1e-6 * _MeV;
  /// @brief Joule expressed in eV
  /// @attention This constant can only be used to convert eV-based energy
  ///            quantities, but not with mechanical units (e.g. kg * m^2 / s^2)
  const double _J = 1. / 1.60217733e-19 * _eV;
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
  const double _T = _kg / (_C * _s);
  const double _N = _kg * _m * _m / (_s * _s);
  /// @}

  /// @name conversion and other constants
  /// @{
  const double _c            = 2.99792458e8 * _m / _s;
  const double _hbar         = 1.05457266e-34 * _J * _s;
  const double _fm_times_GeV = _c * _hbar / (_fm * _GeV);
  const double _kg_per_MeV   = 1. / (_c * _c) * (_kg / _J) / (_kg / _MeV);
  /// @}
}  // namespace units

}  // namespace Acts
#endif  // ACTS_UNITS_HPP
