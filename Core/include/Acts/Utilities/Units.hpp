// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts {

/// @brief Unit and conversion constants
///
/// In order to make sure that always consistent numbers and units are used in
/// calculations, one should make use of the constants defined in this
/// namespace to express the units of numerical values. The following
/// conventions are used:
/// - Input values to Acts must be given as `numerical_value * unit_constant`.
/// - Output values can be converted into the desired unit using
///   `numerical_value / unit_constant`.
/// - Converting between SI and natural units should be done using the template
///   functions `SI2Nat(const double)` and `Nat2SI(const double)`.
///
/// Examples:
/// @code
/// #include "Acts/include/Utilities/Units.hpp"
/// using namespace Acts::units;
/// // specify input variables
/// double momentum = 2.5 * _GeV;
/// double width    = 23 * _cm;
/// double velocity = 345 * _m/_s;
/// double density  = 1.2 * _kg/(_m*_m*_m);
/// double bfield   = 2 * _T;
///
/// // convert output values
/// double x_in_mm   = trackPars.position().x() / _mm;
/// double pt_in_TeV = trackPars.momentum().pT() / _TeV;
/// @endcode
namespace units {

/// @name length units
/// @{
#ifdef DOXYGEN
  constexpr double _m = unspecified;
#else
  constexpr double _m = 1e3;
#endif  // DOXYGEN
  constexpr double _km = 1e3 * _m;
  constexpr double _cm = 1e-2 * _m;
  constexpr double _mm = 1e-3 * _m;
  constexpr double _um = 1e-6 * _m;
  constexpr double _nm = 1e-9 * _m;
  constexpr double _pm = 1e-12 * _m;
  constexpr double _fm = 1e-15 * _m;
  /// Higher orders
  constexpr double _mm2 = _mm * _mm;
/// @}

/// @name mass units
/// @{
#ifdef DOXYGEN
  constexpr double _kg = unspecified;
#else
  constexpr double _kg  = 1e3;
#endif  // DOXYGEN
  constexpr double _g  = 1e-3 * _kg;
  constexpr double _mg = 1e-6 * _kg;
  /// atomic mass unit
  constexpr double _u = 1.660539040e-27 * _kg;
/// @}

/// @name time units
/// @{
#ifdef DOXYGEN
  constexpr double _s = unspecified;
#else
  constexpr double _s   = 1;
#endif  // DOXYGEN
  constexpr double _ms = 1e-3 * _s;
  constexpr double _h  = 3600 * _s;
/// @}

/// @name energy units
/// @{
#ifdef DOXYGEN
  constexpr double _MeV = unspecified;
#else
  constexpr double _MeV = 1e-3;
#endif  // DOXYGEN
  constexpr double _GeV = 1e3 * _MeV;
  constexpr double _TeV = 1e6 * _MeV;
  constexpr double _keV = 1e-3 * _MeV;
  constexpr double _eV  = 1e-6 * _MeV;
/// @}

/// @name charge units
/// @{
#ifdef DOXYGEN
  constexpr double _C = unspecified;
#else
  constexpr double _C   = 1. / 1.60217733e-19;
#endif  // DOXYGEN
  constexpr double _e = 1.60217733e-19 * _C;
  /// Higher orders
  constexpr double _e2 = _e * _e;
  /// @}

  /// @name derived units
  /// @{
  constexpr double _N      = _kg * _m / (_s * _s);
  constexpr double _J      = _N * _m;
  constexpr double _T      = _kg / (_C * _s);
  constexpr double _Gauss  = 1e-4 * _T;
  constexpr double _kGauss = 1e-1 * _T;
  /// @}

  /// @name fundamental physical constants in SI units
  /// @{
  /// speed of light in vacuum
  constexpr double _c = 2.99792458e8 * _m / _s;
  /// reduced Planck constant
  constexpr double _hbar = 1.05457266e-34 * _J * _s;
  /// value of elementary charge in Coulomb
  constexpr double _el_charge = _e / _C;
  /// Higher orders
  constexpr double _c2 = _c * _c;
  constexpr double _c3 = _c * _c * _c;
  constexpr double _c4 = _c2 * _c2;
  /// @}

  /// @cond
  /// @brief internally used conversion constants
  namespace {
    // 1 GeV = 1e9 * e * 1 V = 1.60217733e-10 As * 1 V = 1.60217733e-10 J
    // ==> 1 J = 1 / 1.60217733e-10 GeV
    constexpr double _GeV_per_J = _GeV / (_el_charge * 1e9 * _J);
    // hbar * c = 3.161529298809983e-26 J * m
    // ==> hbar * c * _GeV_per_J = 1.973270523563071e-16 GeV * m
    constexpr double _mm_times_GeV = _c * _hbar * _GeV_per_J;
  }
  /// @endcond

  /// @brief physical quantities for selecting right conversion function
  enum Quantity { MOMENTUM, ENERGY, LENGTH, MASS };

  /// @cond
  template <Quantity>
  double
  SI2Nat(const double);

  template <Quantity>
  double
  Nat2SI(const double);
  /// @endcond

  /// @brief convert energy from SI to natural units
  ///
  /// This function converts the given energy value from SI units to natural
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double E_in_TeV = SI2Nat<ENERGY>(2.3 * _J) / _TeV;
  /// @endcode
  ///
  /// @param[in] E numeric value of energy in SI units
  /// @result numeric value of energy in natural units
  template <>
  double
  SI2Nat<ENERGY>(const double E);

  /// @brief convert energy from natural to SI units
  ///
  /// This function converts the given energy value from natural units to SI
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double E_in_Joule = Nat2SI<ENERGY>(2.3 * _TeV) / _J;
  /// @endcode
  ///
  /// @param[in] E numeric value of energy in natural units
  /// @result numeric value of energy in SI units
  template <>
  double
  Nat2SI<ENERGY>(const double E);

  /// @brief convert length from SI to natural units
  ///
  /// This function converts the given length value from SI units to natural
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double l_per_MeV = SI2Nat<LENGTH>(3 * _um) * _MeV;
  /// @endcode
  ///
  /// @param[in] l numeric value of length in SI units
  /// @result numeric value of length in natural units
  template <>
  double
  SI2Nat<LENGTH>(const double l);

  /// @brief convert length from natural to SI units
  ///
  /// This function converts the given length value from natural units to SI
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double l_in_m = Nat2SI<LENGTH>(1. / (2 * _TeV)) / _m;
  /// @endcode
  ///
  /// @param[in] l numeric value of length in natural units
  /// @result numeric value of length in SI units
  template <>
  double
  Nat2SI<LENGTH>(const double l);

  /// @brief convert momentum from SI to natural units
  ///
  /// This function converts the given momentum value from SI units to natural
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double p_in_GeV = SI2Nat<MOMENTUM>(2 * _N * _s) / _GeV;
  /// @endcode
  ///
  /// @param[in] p numeric value of momentum in SI units
  /// @result numeric value of momentum in natural units
  template <>
  double
  SI2Nat<MOMENTUM>(const double p);

  /// @brief convert momentum from natural to SI units
  ///
  /// This function converts the given momentum value from natural units to SI
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double p_in_Ns = Nat2SI<MOMENTUM>(132 * _GeV) / (_N * _s);
  /// @endcode
  ///
  /// @param[in] p numeric value of momentum in natural units
  /// @result numeric value of momentum in SI units
  template <>
  double
  Nat2SI<MOMENTUM>(const double p);

  /// @brief convert mass from SI to natural units
  ///
  /// This function converts the given mass value from SI units to natural
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double m_in_keV = SI2Nat<MASS>(2 * _g) / _keV;
  /// @endcode
  ///
  /// @param[in] m numeric value of mass in SI units
  /// @result numeric value of mass in natural units
  template <>
  double
  SI2Nat<MASS>(const double m);

  /// @brief convert mass from natural to SI units
  ///
  /// This function converts the given mass value from natural units to SI
  /// units. Example:
  /// @code
  /// #include "Acts/include/Utilities/Units.hpp"
  /// using namespace Acts::units;
  ///
  /// double higgs_in_kg= Nat2SI<MASS>(125 * _GeV) / _kg;
  /// @endcode
  ///
  /// @param[in] m numeric value of mass in natural units
  /// @result numeric value of mass in SI units
  template <>
  double
  Nat2SI<MASS>(const double m);
}  // namespace units

}  // namespace Acts
