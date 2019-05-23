// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

/// Unit definitions and conversions.
///
/// All physical quantities have both a numerical value and a unit. For the
/// computations we always choose a particular unit so we only need to consider
/// the numerical values as such. The chosen base unit for a particular
/// physical quantity, e.g. length, time, or energy, within this code base
/// is called the native unit.
///
/// Here, the following native units are used:
///
/// *   Length is expressed in mm.
/// *   Time is expressed in [speed-of-light * time] == mm. A consequence
///     of this choice is that the speed-of-light expressed in native units
///     is 1.
/// *   Angles are expressed in radian.
/// *   Energy, mass, and momentum are all expressed in GeV (consistent with
///     a speed-of-light == 1).
/// *   Electric charge is expressed in e, i.e. units of the elementary charge.
/// *   The magnetic field is expressed in GeV/(e*mm). The magnetic field
///     connects momentum to length, e.g. in SI units the radius of a charged
///     particle trajectory in a constant magnetic field is given by
///
///         radius = - (momentum / charge) / magnetic-field
///
///     With the chosen magnetic field unit the expression above stays the
///     same and no additional conversion factors are necessary.
///
/// To ensure consistent computations and results the following guidelines
/// **must** be followed when handling physical quantities with units:
///
/// *   All unqualified numerical values, i.e. without a unit, are assumed to
///     be expressed in the relevant native unit, e.g. mm for lengths or GeV
///     for energy/momentum.
/// *   If a variable stores a physical quantity in a specific unit that is
///     not the native unit, clearly mark this in the variable, i.e.
///
///         double momentum = 100.0; // momentum is stored as native unit GeV
///         double momentumInMeV = 10.0; // would be 0.01 in native units
///
/// *   All input values must be given as `numerical_value * unit_constant` or
///     equivalently using the unit literals as `value_unit`. The resulting
///     unqualified numerical value will be automatically converted to the
///     native unit.
/// *   To output an unqualified numerical value in the native units as a
///     numerical value in a specific unit divide by the unit constants as
///     `numerical_value / unit_constant` or using the unit literals as
///     `value / 1_unit`.
///
/// Examples:
///
///     #include "Acts/include/Utilities/Units.hpp"
///     using namespace Acts::UnitLiterals;
///
///     // define input values w/ units (via unit constants)
///     double width    = 12 * Acts::UnitConstants::mm;
///     double mmuon    = 105.7 * Acts::UnitConstants::MeV;
///     // define input values w/ units (via unit user literals)
///     double length   = 23_cm;
///     double time     = 1214.2_ns;
///     double angle    = 123_degree;
///     double momentum = 2.5_TeV;
///     double mass     = 511_keV;
///     double velocity = 345_m / 1_s;
///     double bfield   = 3.9_T;
///
///     // convert output values (via unit constants)
///     doube t_in_ns    = trackPars.time() / Acts::UnitConstants::ns;
///     // convert output values (via unit user literals)
///     double x_in_mm   = trackPars.position().x() / 1_mm;
///     double pt_in_TeV = trackPars.momentum().pT() / 1_TeV;
///

namespace UnitConstants {
// Length, native unit mm
constexpr double fm = 1e-12;
constexpr double pm = 1e-9;
constexpr double um = 1e-3;
constexpr double nm = 1e-6;
constexpr double mm = 1.0;
constexpr double cm = 10.0;
constexpr double m = 1e3;
constexpr double km = 1e6;
// Time, native unit mm = [speed-of-light * time] = mm/s * s
constexpr double s = 299792458000.0;
constexpr double fs = 1e-15 * s;
constexpr double ps = 1e-12 * s;
constexpr double ns = 1e-9 * s;
constexpr double us = 1e-6 * s;
constexpr double ms = 1e-3 * s;
constexpr double min = 60.0 * s;
constexpr double h = 3600.0 * s;
// Angles, native unit radian
constexpr double mrad = 1e-3;
constexpr double rad = 1.0;
constexpr double degree = 0.017453292519943295;  // pi / 180
// Energy/mass/momentum, native unit GeV
constexpr double eV = 1e-9;
constexpr double keV = 1e-6;
constexpr double MeV = 1e-3;
constexpr double GeV = 1.0;
constexpr double TeV = 1e3;
//     1eV/c² == 1.782662e-36kg
//    1GeV/c² == 1.782662e-27kg
// ->     1kg == (1/1.782662e-27)GeV/c²
// ->      1g == (1/(1e3*1.782662e-27))GeV/c²
constexpr double g = 1.0 / 1.782662e-24;
constexpr double kg = 1.0 / 1.782662e-27;
// Charge, native unit e (elementary charge)
constexpr double e = 1.0;
constexpr double C = 1.602176634e19;
// Magnetic field, native unit GeV/(e*mm)
constexpr double T = 0.000299792458;  // equivalent to c in appropriate SI units
constexpr double Gauss = 1e-4 * T;
constexpr double kGauss = 1e-1 * T;
}  // namespace UnitConstants

namespace UnitLiterals {
// define user literal functions for the given unit constant
#define ACTS_DEFINE_UNIT_LITERAL(name)                        \
  constexpr double operator"" _##name(long double x) {        \
    return ::Acts::UnitConstants::name * x;                   \
  }                                                           \
  constexpr double operator"" _##name(unsigned long long x) { \
    return ::Acts::UnitConstants::name * x;                   \
  }
ACTS_DEFINE_UNIT_LITERAL(fm)
ACTS_DEFINE_UNIT_LITERAL(pm)
ACTS_DEFINE_UNIT_LITERAL(nm)
ACTS_DEFINE_UNIT_LITERAL(um)
ACTS_DEFINE_UNIT_LITERAL(mm)
ACTS_DEFINE_UNIT_LITERAL(cm)
ACTS_DEFINE_UNIT_LITERAL(m)
ACTS_DEFINE_UNIT_LITERAL(km)
ACTS_DEFINE_UNIT_LITERAL(fs)
ACTS_DEFINE_UNIT_LITERAL(ps)
ACTS_DEFINE_UNIT_LITERAL(ns)
ACTS_DEFINE_UNIT_LITERAL(us)
ACTS_DEFINE_UNIT_LITERAL(ms)
ACTS_DEFINE_UNIT_LITERAL(s)
ACTS_DEFINE_UNIT_LITERAL(min)
ACTS_DEFINE_UNIT_LITERAL(h)
ACTS_DEFINE_UNIT_LITERAL(mrad)
ACTS_DEFINE_UNIT_LITERAL(rad)
ACTS_DEFINE_UNIT_LITERAL(degree)
ACTS_DEFINE_UNIT_LITERAL(eV)
ACTS_DEFINE_UNIT_LITERAL(keV)
ACTS_DEFINE_UNIT_LITERAL(MeV)
ACTS_DEFINE_UNIT_LITERAL(GeV)
ACTS_DEFINE_UNIT_LITERAL(TeV)
ACTS_DEFINE_UNIT_LITERAL(g)
ACTS_DEFINE_UNIT_LITERAL(kg)
ACTS_DEFINE_UNIT_LITERAL(e)
ACTS_DEFINE_UNIT_LITERAL(C)
ACTS_DEFINE_UNIT_LITERAL(T)
ACTS_DEFINE_UNIT_LITERAL(Gauss)
ACTS_DEFINE_UNIT_LITERAL(kGauss)
// not needed anymore. undef to prevent littering the namespace
#undef ACTS_DEFINE_UNIT_LITERAL
}  // namespace UnitLiterals

/// Legacy namespace for backward-compatibility
namespace units {

/// @name length units
/// @{
constexpr double _km = UnitConstants::km;
constexpr double _m = UnitConstants::m;
constexpr double _cm = UnitConstants::cm;
constexpr double _mm = UnitConstants::mm;
constexpr double _um = UnitConstants::um;
constexpr double _nm = UnitConstants::nm;
constexpr double _pm = UnitConstants::pm;
constexpr double _fm = UnitConstants::fm;
/// Higher orders
constexpr double _mm2 = _mm * _mm;
/// @}

/// @name time units
/// @{
constexpr double _h = UnitConstants::h;
constexpr double _s = UnitConstants::s;
constexpr double _ms = UnitConstants::ms;
/// @}

/// @name energy/momentum/mass units
/// @{
constexpr double _GeV = UnitConstants::GeV;
constexpr double _TeV = UnitConstants::TeV;
constexpr double _MeV = UnitConstants::MeV;
constexpr double _keV = UnitConstants::keV;
constexpr double _eV = UnitConstants::eV;
constexpr double _g = UnitConstants::g;
constexpr double _kg = UnitConstants::kg;
/// @}

/// @name charge units
/// @{
constexpr double _C = UnitConstants::C;
constexpr double _e = UnitConstants::e;
/// Higher orders
constexpr double _e2 = _e * _e;
/// @}

/// @name magnetic field units
/// @{
constexpr double _T = UnitConstants::T;
constexpr double _Gauss = UnitConstants::Gauss;
constexpr double _kGauss = UnitConstants::kGauss;
/// @}

/// @brief physical quantities for selecting right conversion function
enum Quantity { MOMENTUM, ENERGY, LENGTH, MASS };

/// @cond
template <Quantity>
double SI2Nat(const double);

template <Quantity>
double Nat2SI(const double);
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
double SI2Nat<ENERGY>(const double E);

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
double Nat2SI<ENERGY>(const double E);

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
double SI2Nat<LENGTH>(const double l);

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
double Nat2SI<LENGTH>(const double l);

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
double SI2Nat<MOMENTUM>(const double p);

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
double Nat2SI<MOMENTUM>(const double p);

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
double SI2Nat<MASS>(const double m);

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
double Nat2SI<MASS>(const double m);
}  // namespace units

}  // namespace Acts
