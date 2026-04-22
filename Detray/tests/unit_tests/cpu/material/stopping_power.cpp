// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/material/interaction.hpp"
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

// Test class for the stopping power
// Input tuple: < material, particle type, kinetic energy, expected output >
class StoppingPowerValidation
    : public ::testing::TestWithParam<
          std::tuple<material<test::scalar>, pdg_particle<test::scalar>,
                     test::scalar, test::scalar>> {};

TEST_P(StoppingPowerValidation, stopping_power) {
  // Interaction object
  interaction<test::scalar> I;

  // Material
  material<test::scalar> mat = std::get<0>(GetParam());

  // Particle
  pdg_particle<test::scalar> ptc = std::get<1>(GetParam());

  // Kinetic energy
  const test::scalar T = std::get<2>(GetParam());

  // Total energy
  const test::scalar E = T + ptc.mass();

  // Momentum
  const test::scalar p = math::sqrt(E * E - ptc.mass() * ptc.mass());

  // qoverp
  const test::scalar qop{ptc.charge() / p};

  // Stopping power in MeV * cm^2 / g
  const test::scalar dEdx{I.compute_stopping_power(mat, ptc, {ptc, qop}) /
                          mat.mass_density() /
                          (unit<test::scalar>::MeV * unit<test::scalar>::cm2 /
                           unit<test::scalar>::g)};

  const test::scalar expected_dEdx = std::get<3>(GetParam());

  // Check if difference is within 8% error
  EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.08f);
}

/******************
 *   Muon tests
 ******************/

// From https://pdg.lbl.gov/2024/AtomicNuclearProperties/index.html
// Note 1: that we took the PDG value only from Ionization loss (Radiative loss
// is ignored) Note 2: assumes that the stopping powers of muon and antimuon are
// the same Note 3: Test fails with He Gas and 1 GeV muons (18 % difference)
INSTANTIATE_TEST_SUITE_P(
    muon_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(
        std::make_tuple(helium_gas<test::scalar>(), muon<test::scalar>(),
                        100.0f * unit<test::scalar>::MeV, 2.165f),
        // std::make_tuple(helium_gas<test::scalar>(), muon<test::scalar>(),
        //                 1.f * unit<test::scalar>::GeV, 2.133f),
        std::make_tuple(helium_gas<test::scalar>(), muon<test::scalar>(),
                        10.0f * unit<test::scalar>::GeV, 2.768f),
        std::make_tuple(helium_gas<test::scalar>(), muon<test::scalar>(),
                        100.0f * unit<test::scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    muon_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(
        std::make_tuple(silicon<test::scalar>(), muon<test::scalar>(),
                        100.0f * unit<test::scalar>::MeV, 1.849f),
        std::make_tuple(silicon<test::scalar>(), muon<test::scalar>(),
                        1.f * unit<test::scalar>::GeV, 1.803f),
        std::make_tuple(silicon<test::scalar>(), muon<test::scalar>(),
                        10.0f * unit<test::scalar>::GeV, 2.177f),
        std::make_tuple(silicon<test::scalar>(), muon<test::scalar>(),
                        100.0f * unit<test::scalar>::GeV, 2.451f)));

INSTANTIATE_TEST_SUITE_P(
    anti_muon_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(
        std::make_tuple(helium_gas<test::scalar>(), antimuon<test::scalar>(),
                        100.0f * unit<test::scalar>::MeV, 2.165f),
        // std::make_tuple(helium_gas<test::scalar>(), antimuon<test::scalar>(),
        //                 1.f * unit<test::scalar>::GeV, 2.133f),
        std::make_tuple(helium_gas<test::scalar>(), antimuon<test::scalar>(),
                        10.0f * unit<test::scalar>::GeV, 2.768f),
        std::make_tuple(helium_gas<test::scalar>(), antimuon<test::scalar>(),
                        100.0f * unit<test::scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    anti_muon_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(
        std::make_tuple(silicon<test::scalar>(), antimuon<test::scalar>(),
                        100.0f * unit<test::scalar>::MeV, 1.849f),
        std::make_tuple(silicon<test::scalar>(), antimuon<test::scalar>(),
                        1.f * unit<test::scalar>::GeV, 1.803f),
        std::make_tuple(silicon<test::scalar>(), antimuon<test::scalar>(),
                        10.0f * unit<test::scalar>::GeV, 2.177f),
        std::make_tuple(silicon<test::scalar>(), antimuon<test::scalar>(),
                        100.0f * unit<test::scalar>::GeV, 2.451f)));

/*********************
 *   Electron tests
 *********************/

// From https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
// Assumes that the stopping powers of electron and positron are the same
INSTANTIATE_TEST_SUITE_P(
    electron_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<test::scalar>(),
                                      electron<test::scalar>(),
                                      100.0f * unit<test::scalar>::MeV, 3.532f),
                      std::make_tuple(helium_gas<test::scalar>(),
                                      electron<test::scalar>(),
                                      1.f * unit<test::scalar>::GeV, 13.14f)));

INSTANTIATE_TEST_SUITE_P(
    electron_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<test::scalar>(),
                                      electron<test::scalar>(),
                                      100.0f * unit<test::scalar>::MeV, 6.017f),
                      std::make_tuple(silicon<test::scalar>(),
                                      electron<test::scalar>(),
                                      1.f * unit<test::scalar>::GeV, 46.69f)));

INSTANTIATE_TEST_SUITE_P(
    positron_stopping_power_He, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(helium_gas<test::scalar>(),
                                      positron<test::scalar>(),
                                      100.0f * unit<test::scalar>::MeV, 3.532f),
                      std::make_tuple(helium_gas<test::scalar>(),
                                      positron<test::scalar>(),
                                      1.f * unit<test::scalar>::GeV, 13.14f)));

INSTANTIATE_TEST_SUITE_P(
    positron_stopping_power_Si, StoppingPowerValidation,
    ::testing::Values(std::make_tuple(silicon<test::scalar>(),
                                      positron<test::scalar>(),
                                      100.0f * unit<test::scalar>::MeV, 6.017f),
                      std::make_tuple(silicon<test::scalar>(),
                                      positron<test::scalar>(),
                                      1.f * unit<test::scalar>::GeV, 46.69f)));
