// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

using namespace Acts::UnitLiterals;

BOOST_AUTO_TEST_SUITE(unit_conversions)

BOOST_AUTO_TEST_CASE(length_conversions) {
  CHECK_CLOSE_REL(1_m, 1e-3_km, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e2_cm, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e3_mm, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e6_um, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e9_nm, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e12_pm, 1e-15);
  CHECK_CLOSE_REL(1_m, 1e15_fm, 1e-15);
}

BOOST_AUTO_TEST_CASE(area_conversions) {
  CHECK_CLOSE_REL(1_mm * 1_mm, 1_mm2, 1e-15);
  CHECK_CLOSE_REL(1_mm * 1_mm, 0.01_cm2, 1e-15);
  CHECK_CLOSE_REL(1_mm * 1_mm, 0.000001_m2, 1e-15);
  CHECK_CLOSE_REL(1_cm * 1_cm, 100_mm2, 1e-15);
  CHECK_CLOSE_REL(1_cm * 1_cm, 1_cm2, 1e-15);
  CHECK_CLOSE_REL(1_cm * 1_cm, 0.0001_m2, 1e-15);
  CHECK_CLOSE_REL(1_m * 1_m, 1000000_mm2, 1e-15);
  CHECK_CLOSE_REL(1_m * 1_m, 10000_cm2, 1e-15);
  CHECK_CLOSE_REL(1_m * 1_m, 1_m2, 1e-15);
}

BOOST_AUTO_TEST_CASE(volume_conversions) {
  CHECK_CLOSE_REL(1_mm * 1_mm * 1_mm, 1_mm3, 1e-15);
  CHECK_CLOSE_REL(1_mm2 * 1_mm, 1_mm3, 1e-15);
  CHECK_CLOSE_REL(1_cm * 1_cm * 1_cm, 1_cm3, 1e-15);
  CHECK_CLOSE_REL(1_cm2 * 1_cm, 1_cm3, 1e-15);
  CHECK_CLOSE_REL(1_m * 1_m * 1_m, 1_m3, 1e-15);
  CHECK_CLOSE_REL(1_m2 * 1_m, 1_m3, 1e-15);
}

BOOST_AUTO_TEST_CASE(time_conversions) {
  CHECK_CLOSE_REL(1_h, 60_min, 1e-15);
  CHECK_CLOSE_REL(1_h, 3600_s, 1e-15);
  CHECK_CLOSE_REL(1_min, 60_s, 1e-15);
  CHECK_CLOSE_REL(1_s, 1e3_ms, 1e-15);
  CHECK_CLOSE_REL(1_s, 1e6_us, 1e-15);
  CHECK_CLOSE_REL(1_s, 1e9_ns, 1e-15);
  CHECK_CLOSE_REL(1_s, 1e12_ps, 1e-15);
  CHECK_CLOSE_REL(1_s, 1e15_fs, 1e-15);
}

BOOST_AUTO_TEST_CASE(angle_conversions) {
  CHECK_CLOSE_REL(45_degree, M_PI / 4 * 1_rad, 1e-15);
  CHECK_CLOSE_REL(90_degree, M_PI / 2 * 1_rad, 1e-15);
  CHECK_CLOSE_REL(180_degree, M_PI * 1_rad, 1e-15);
  CHECK_CLOSE_REL(360_degree, 2 * M_PI * 1_rad, 1e-15);
  CHECK_CLOSE_REL(1_mm / 1_m, 1_mrad, 1e-15);
  CHECK_CLOSE_REL(1_um / 1_mm, 1_mrad, 1e-15);
}

BOOST_AUTO_TEST_CASE(energy_conversions) {
  CHECK_CLOSE_REL(1_MeV, 1e6_eV, 1e-15);
  CHECK_CLOSE_REL(1_MeV, 1e3_keV, 1e-15);
  CHECK_CLOSE_REL(1_MeV, 1e-3_GeV, 1e-15);
  CHECK_CLOSE_REL(1_MeV, 1e-6_TeV, 1e-15);
}

BOOST_AUTO_TEST_CASE(mass_conversions) {
  CHECK_CLOSE_REL(1_kg, 1000_g, 1e-15);
  CHECK_CLOSE_REL(0.001_kg, 1_g, 1e-15);
  CHECK_CLOSE_REL(1_u, 931.49410242_MeV, 1e-15);
  CHECK_CLOSE_REL(1_u, 1.66053906660e-24_g, 1e-7);
}

BOOST_AUTO_TEST_CASE(mass_energy_conversions) {
  // always assumes c == 1
  CHECK_CLOSE_REL(1.782662e-36_kg, 1_eV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-33_kg, 1_keV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-30_kg, 1_MeV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-27_kg, 1_GeV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-33_g, 1_eV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-30_g, 1_keV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-27_g, 1_MeV, 1e-15);
  CHECK_CLOSE_REL(1.782662e-24_g, 1_GeV, 1e-15);
}

BOOST_AUTO_TEST_CASE(charge_conversions) {
  CHECK_CLOSE_REL(1_C, 1.602176634e19_e, 1e-15);
}

BOOST_AUTO_TEST_CASE(magnetic_field_conversions) {
  CHECK_CLOSE_REL(10_kGauss, 1_T, 1e-15);
  CHECK_CLOSE_REL(1_kGauss, 1000_Gauss, 1e-15);
}

BOOST_AUTO_TEST_CASE(momentum_to_radius_conversion) {
  // note: no conversion factors necessary
  CHECK_CLOSE_REL(1_GeV / (1_e * 1_T), 3.3336_m, 1e-3);
  CHECK_CLOSE_REL(1_GeV / (1_e * 2_T), 166.8_cm, 1e-3);
  CHECK_CLOSE_REL(1_GeV / (2_e * 1_T), 166.8_cm, 1e-3);
  CHECK_CLOSE_REL(1_GeV / (1_e * 4_T), 83.39_cm, 1e-3);
  CHECK_CLOSE_REL(1_GeV / (2_e * 2_T), 83.39_cm, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
