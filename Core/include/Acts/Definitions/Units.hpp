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

/// @verbatim embed:rst:leading-slashes
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
///   .. math::
///
///      radius = - (momentum / charge) / field
///
///   With the chosen magnetic field unit the expression above stays the
///   same and no additional conversion factors are necessary.
/// - Amount of substance is expressed in mol.
///
/// Depending on the context a physical quantity might not be given in the
/// native units. In this case we need to convert to the native unit first
/// before the value can be used. The necessary conversion factors are defined
/// in  ``Acts::UnitConstants``. Multiplying a value with the unit constant
/// produces the equivalent value in the native unit, e.g.
///
/// .. code-block:: cpp
///
///    long double length = 1 * Acts::UnitConstants::m;       // length ==
///    1000.0 long double momentum = 100 * Acts::UnitConstants::MeV; // momentum
///    == 0.1
///
/// The conversion of a value in native units into the equivalent value in a
/// specific other unit is computed by dividing with the relevant unit, e.g.
///
/// .. code-block:: cpp
///
///    long double length = 10.0;                               // native units
///    mm long double lengthInM = length / Acts::UnitConstants::m; // == 0.01;
///
/// To further simplify the usage, physical quantities can also be expressed via
/// `C++ user literals
/// <https://en.cppreference.com/w/cpp/language/user_literal>`_ defined in
/// ``Acts::UnitLiterals``. This allows us to express quantities in a concise
/// way:
///
/// .. code-block:: cpp
///
///    using namespace Acts::UnitLiterals;
///
///    long double length = 1_m;                     // == 1000.0
///    long double momentum = 1.25_TeV;              // == 1250.0
///    long double lengthInUm = length / 1_um;       // == 1000000.0
///    long double momentumInMeV = momentum / 1_MeV; // == 1250000.0
///
/// .. warning::
///    Since using user-defined literals requires a namespace import of
///    ``Acts::UnitLiterals`` it should not be used in headers since it would
///    (accidentally) modify the namespace wherever the header is included.
///
/// To ensure consistent computations and results the following guidelines
/// **must** be followed when handling physical quantities with units:
///
/// - All unqualified numerical values, i.e. without a unit, are assumed to
///   be expressed in the relevant native unit, e.g. mm for lengths or GeV
///   for energy/momentum.
/// - If a variable stores a physical quantity in a specific unit that is
///   not the native unit, clearly mark this in the variable, i.e.
///
///   .. code-block:: cpp
///
///      long double momentum = 100.0; // momentum is stored as native unit GeV
///      long double momentumInMeV = 10.0; // would be 0.01 in native units
///
/// - All input values must be given as ``numerical_value * unit_constant`` or
///   equivalently using the unit literals as ``value_unit``. The resulting
///   unqualified numerical value will be automatically converted to the
///   native unit.
/// - To output an unqualified numerical value in the native units as a
///   numerical value in a specific unit divide by the unit constants as
///   ``numerical_value / unit_constant`` or using the unit literals as
///   ``value / 1_unit``.
///
/// Examples:
///
/// .. code-block:: cpp
///
///    #include <Acts/include/Definitions/Units.hpp>
///    using namespace Acts::UnitLiterals;
///
///    // define input values w/ units (via unit constants)
///    long double width    = 12 * Acts::UnitConstants::mm;
///    long double mmuon    = 105.7 * Acts::UnitConstants::MeV;
///    // define input values w/ units (via unit user literals)
///    long double length   = 23_cm;
///    long double time     = 1214.2_ns;
///    long double angle    = 123_degree;
///    long double momentum = 2.5_TeV;
///    long double mass     = 511_keV;
///    long double velocity = 345_m / 1_s;
///    long double bfield   = 3.9_T;
///
///    // convert output values (via unit constants)
///    long double t_in_ns    = trackPars.time() / Acts::UnitConstants::ns;
///    // convert output values (via unit user literals)
///    long double x_in_mm   = trackPars.position()[ePos0] / 1_mm;
///    long double p_in_TeV = trackPars.absoluteMomentum() / 1_TeV;
///
/// @endverbatim

/// @note A helper script is available in
///   `Core/scripts/print_units_physical_constants.py` to validate some of the
///   numerical values.

namespace UnitConstants {
// Length, native unit mm
constexpr long double mm = 1.0;
constexpr long double fm = 1e-12 * mm;
constexpr long double pm = 1e-9 * mm;
constexpr long double nm = 1e-6 * mm;
constexpr long double um = 1e-3 * mm;
constexpr long double cm = 1e1 * mm;
constexpr long double m = 1e3 * mm;
constexpr long double km = 1e6 * mm;
// Shortcuts for commonly used area and volume units. This intentionally
// contains not all possible combinations to avoid cluttering the namespace.
// Missing area or volume units can always be defined on the fly using the
// existing length units e.g. 1fm³ -> 1fm * 1fm * 1fm
// Area, native unit mm²
constexpr long double mm2 = mm * mm;
constexpr long double cm2 = cm * cm;
constexpr long double m2 = m * m;
// Volume, native unit mm³
constexpr long double mm3 = mm * mm * mm;
constexpr long double cm3 = cm * cm * cm;
constexpr long double m3 = m * m * m;
// Time, native unit mm = [speed-of-light * time] = mm/s * s
/// @note Depends on speed of light in SI units
constexpr long double s = 299792458000.0;  // = 299792458.0 * (m / 1.0) * 1.0;
constexpr long double fs = 1e-15 * s;
constexpr long double ps = 1e-12 * s;
constexpr long double ns = 1e-9 * s;
constexpr long double us = 1e-6 * s;
constexpr long double ms = 1e-3 * s;
constexpr long double min = 60.0 * s;
constexpr long double h = 3600.0 * s;
// Angles, native unit radian
constexpr long double mrad = 1e-3;
constexpr long double rad = 1.0;
constexpr long double degree = std::numbers::pi / 180. / rad;
// Energy/mass/momentum, native unit GeV
constexpr long double GeV = 1.0;
constexpr long double eV = 1e-9 * GeV;
constexpr long double keV = 1e-6 * GeV;
constexpr long double MeV = 1e-3 * GeV;
constexpr long double TeV = 1e3 * GeV;
constexpr long double J = 6241509074.460763 * GeV;
/// atomic mass unit u
constexpr long double u = 0.93149410242;
//     1eV/c² == 1.782662e-36kg
//    1GeV/c² == 1.782662e-27kg
// ->     1kg == (1/1.782662e-27)GeV/c²
// ->      1g == (1/(1e3*1.782662e-27))GeV/c²
constexpr long double g = 1.0 / 1.782662e-24;
constexpr long double kg = 1.0 / 1.782662e-27;
/// Charge, native unit e (elementary charge)
constexpr long double e = 1.0;
/// Magnetic field, native unit (eV*s)/(e*m²)
/// @note Depends on speed of light in SI units
constexpr long double T = 0.000299792458;  // = eV * s / (e * m2);
constexpr long double Gauss = 1e-4 * T;
constexpr long double kGauss = 1e-1 * T;
/// Amount of substance, native unit mol
constexpr long double mol = 1.0;
}  // namespace UnitConstants

namespace UnitLiterals {
// define user literal functions for the given unit constant
#define ACTS_DEFINE_UNIT_LITERAL(name)                            \
  constexpr long double operator""_##name(long long double x) {   \
    return ::Acts::UnitConstants::name * x;                       \
  }                                                               \
  constexpr long double operator""_##name(unsigned long long x) { \
    return ::Acts::UnitConstants::name * x;                       \
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
constexpr long double c = 1.0;
/// Reduced Planck constant h/2*pi.
///
/// Computed from CODATA 2018 constants to long double precision.
constexpr long double hbar =
    6.582119569509066e-25 * UnitConstants::GeV * UnitConstants::s;
}  // namespace PhysicalConstants

}  // namespace Acts
