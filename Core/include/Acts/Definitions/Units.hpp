// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <numbers>

namespace Acts {

/// @namespace Acts::UnitConstants
/// @brief Constants and helper literals for physical units.
/// All physical quantities have both a numerical value and a unit. For the
/// computations we always choose a particular unit for a given physical
/// quantity so we only need to consider the numerical values as such. The
/// chosen base unit for a particular physical dimension, e.g. length, time, or
/// energy, within this code base is called the native unit. The base units
/// should be chosen such that they are internally consistent, require the least
/// amount of explicit conversion factors (ideally none at all), and have
/// typical numerical values close to unity to reduce numerical issues.
///
/// Here, the following native units are used:
///
/// - Length is expressed in mm.
/// - Time is expressed in [speed-of-light * time] == mm. A consequence
///   of this choice is that the speed-of-light expressed in native units
///   is 1.
/// - Angles are expressed in radian.
/// - Energy, mass, and momentum are all expressed in GeV (consistent with
///   speed-of-light == 1).
/// - Electric charge is expressed in e, i.e. units of the elementary charge.
/// - The magnetic field is expressed in GeV/(e*mm). The magnetic field
///   connects momentum to length, e.g. in SI units the radius of a charged
///   particle trajectory in a constant magnetic field is given by
///
///   @f$
///      \text{radius} = - \frac{\text{momentum} }{ \text{charge} } / \text{field}
///   @f$
///
/// With the chosen magnetic field unit the expression above stays the
/// same and no additional conversion factors are necessary.
///
/// Depending on the context a physical quantity might not be given in the
/// native units. In this case we need to convert to the native unit first
/// before the value can be used. The necessary conversion factors are defined
/// in  @ref Acts::UnitConstants. Multiplying a value with the unit constant
/// produces the equivalent value in the native unit, e.g.
///
/// @snippet{trimleft} examples/units.cpp Using Unit Constants
///
/// To further simplify the usage, physical quantities can also be expressed
/// via [C++ user
/// literals](https://en.cppreference.com/w/cpp/language/user_literal) defined
/// in @ref Acts::UnitLiterals. This allows us to express quantities in a
/// concise way:
///
/// @snippet{trimleft} examples/units.cpp Using Unit Literals
///
/// @warning Since using user-defined literals requires a namespace import of
///          @ref Acts::UnitLiterals it should not be used in headers since it
///          would (accidentally) modify the namespace wherever the header is
///          included.
///
/// To ensure consistent computations and results the following guidelines
/// **must** be followed when handling physical quantities with units:
///
/// @snippet{trimleft} examples/units.cpp Unit Best Practices
///
/// Here's a comprehensive example showing various ways to work with units:
///
/// @snippet{trimleft} examples/units.cpp Comprehensive Units Example
///
/// Converting output values from native units:
///
/// @snippet{trimleft} examples/units.cpp Converting Output Values
///
/// @note A helper script is available in
///   `Core/scripts/print_units_physical_constants.py` to validate some of the
///   numerical values.

namespace UnitConstants {
// Length, native unit mm
/// Millimeter - native unit for length
constexpr double mm = 1.0;
/// Femtometer - 1e-15 meter
constexpr double fm = 1e-12 * mm;
/// Picometer - 1e-12 meter
constexpr double pm = 1e-9 * mm;
/// Nanometer - 1e-9 meter
constexpr double nm = 1e-6 * mm;
/// Micrometer - 1e-6 meter
constexpr double um = 1e-3 * mm;
/// Centimeter - 1e-2 meter
constexpr double cm = 1e1 * mm;
/// Meter
constexpr double m = 1e3 * mm;
/// Kilometer - 1e3 meter
constexpr double km = 1e6 * mm;
// Shortcuts for commonly used area and volume units. This intentionally
// contains not all possible combinations to avoid cluttering the namespace.
// Missing area or volume units can always be defined on the fly using the
// existing length units e.g. 1fm³ -> 1fm * 1fm * 1fm
// Area, native unit mm²
/// Square millimeter - native unit for area
constexpr double mm2 = mm * mm;
/// Square centimeter
constexpr double cm2 = cm * cm;
/// Square meter
constexpr double m2 = m * m;
// Volume, native unit mm³
/// Cubic millimeter - native unit for volume
constexpr double mm3 = mm * mm * mm;
/// Cubic centimeter
constexpr double cm3 = cm * cm * cm;
/// Cubic meter
constexpr double m3 = m * m * m;
// Time, native unit mm = [speed-of-light * time] = mm/s * s
/// @note Depends on speed of light in SI units
constexpr double s = 299792458000.0;  // = 299792458.0 * (m / 1.0) * 1.0;
/// Femtosecond - 1e-15 second
constexpr double fs = 1e-15 * s;
/// Picosecond - 1e-12 second
constexpr double ps = 1e-12 * s;
/// Nanosecond - 1e-9 second
constexpr double ns = 1e-9 * s;
/// Microsecond - 1e-6 second
constexpr double us = 1e-6 * s;
/// Millisecond - 1e-3 second
constexpr double ms = 1e-3 * s;
/// Minute - 60 seconds
constexpr double min = 60.0 * s;
/// Hour - 3600 seconds
constexpr double h = 3600.0 * s;
// Angles, native unit radian
/// Milliradian - 1e-3 radian
constexpr double mrad = 1e-3;
/// Radian - native unit for angle
constexpr double rad = 1.0;
/// Degree - pi/180 radians
constexpr double degree = std::numbers::pi / 180. / rad;
// Energy/mass/momentum, native unit GeV
/// Gigaelectronvolt - native unit for energy/mass/momentum
constexpr double GeV = 1.0;
/// Electronvolt - 1e-9 GeV
constexpr double eV = 1e-9 * GeV;
/// Kiloelectronvolt - 1e-6 GeV
constexpr double keV = 1e-6 * GeV;
/// Megaelectronvolt - 1e-3 GeV
constexpr double MeV = 1e-3 * GeV;
/// Teraelectronvolt - 1e3 GeV
constexpr double TeV = 1e3 * GeV;
/// Joule in GeV
constexpr double J = 6241509074.460763 * GeV;
/// atomic mass unit u
constexpr double u = 0.93149410242;
/// Gram in GeV/c²
/// @note 1eV/c² == 1.782662e-36kg
///      1GeV/c² == 1.782662e-27kg
///   ->     1kg == (1/1.782662e-27)GeV/c²
///   ->      1g == (1/(1e3*1.782662e-27))GeV/c²
constexpr double g = 1.0 / 1.782662e-24;
/// Kilogram in GeV/c²
constexpr double kg = 1.0 / 1.782662e-27;
/// Charge, native unit e (elementary charge)
constexpr double e = 1.0;
/// Magnetic field, native unit (eV*s)/(e*m²)
/// @note Depends on speed of light in SI units
constexpr double T = 0.000299792458;  // = eV * s / (e * m2);
/// Gauss - 1e-4 Tesla
constexpr double Gauss = 1e-4 * T;
/// Kilogauss - 1e-1 Tesla
constexpr double kGauss = 1e-1 * T;
/// Amount of substance, native unit mol
constexpr double mol = 1.0;
}  // namespace UnitConstants

/// @brief Namespace for user-defined literals for physical units. See @ref
///        UnitConstants for details.
namespace UnitLiterals {
// define user literal functions for the given unit constant
#define ACTS_DEFINE_UNIT_LITERAL(name)                       \
  constexpr double operator""_##name(long double x) {        \
    return ::Acts::UnitConstants::name * x;                  \
  }                                                          \
  constexpr double operator""_##name(unsigned long long x) { \
    return ::Acts::UnitConstants::name * x;                  \
  }
ACTS_DEFINE_UNIT_LITERAL(fm)
ACTS_DEFINE_UNIT_LITERAL(pm)
ACTS_DEFINE_UNIT_LITERAL(nm)
ACTS_DEFINE_UNIT_LITERAL(um)
ACTS_DEFINE_UNIT_LITERAL(mm)
ACTS_DEFINE_UNIT_LITERAL(cm)
ACTS_DEFINE_UNIT_LITERAL(m)
ACTS_DEFINE_UNIT_LITERAL(km)
ACTS_DEFINE_UNIT_LITERAL(mm2)
ACTS_DEFINE_UNIT_LITERAL(cm2)
ACTS_DEFINE_UNIT_LITERAL(m2)
ACTS_DEFINE_UNIT_LITERAL(mm3)
ACTS_DEFINE_UNIT_LITERAL(cm3)
ACTS_DEFINE_UNIT_LITERAL(m3)
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
ACTS_DEFINE_UNIT_LITERAL(J)
ACTS_DEFINE_UNIT_LITERAL(u)
ACTS_DEFINE_UNIT_LITERAL(g)
ACTS_DEFINE_UNIT_LITERAL(kg)
ACTS_DEFINE_UNIT_LITERAL(e)
ACTS_DEFINE_UNIT_LITERAL(T)
ACTS_DEFINE_UNIT_LITERAL(Gauss)
ACTS_DEFINE_UNIT_LITERAL(kGauss)
ACTS_DEFINE_UNIT_LITERAL(mol)
// not needed anymore. undef to prevent littering the namespace
#undef ACTS_DEFINE_UNIT_LITERAL
}  // namespace UnitLiterals

/// Physical constants in native units.
///
/// Unit constants are intentionally not listed.
namespace PhysicalConstants {
// Speed of light
/// Speed of light in vacuum - native unit (dimensionless)
constexpr double c = 1.0;
/// Reduced Planck constant h/2*pi.
///
/// Computed from CODATA 2018 constants to double precision.
constexpr double hbar =
    6.582119569509066e-25 * UnitConstants::GeV * UnitConstants::s;
}  // namespace PhysicalConstants

}  // namespace Acts
