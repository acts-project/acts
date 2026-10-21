// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/definitions/pdg_particle.hpp"
#include "detray/material/interaction.hpp"
#include "detray/material/material.hpp"
#include "detray/material/predefined_materials.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using scalar = test::scalar;

// Test class for energy loss with Bremsstrahlung
// Input tuple: < material, particle type, kinetic energy, expected output >
class EnergyLossBremsValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, pdg_particle<scalar>, scalar, scalar>> {
};

// This tests the material functionalities
TEST_P(EnergyLossBremsValidation, bremsstrahlung) {
  // Interaction object
  interaction<scalar> I;

  // Material
  material<scalar> mat = std::get<0>(GetParam());

  // Particle
  pdg_particle<scalar> ptc = std::get<1>(GetParam());

  // Kinetic energy
  const scalar T = std::get<2>(GetParam());

  // Total energy
  const scalar E = T + ptc.mass();

  // Momentum
  const scalar p = math::sqrt(E * E - ptc.mass() * ptc.mass());

  // qoverp
  const scalar qop{ptc.charge() / p};

  // Bremsstrahlung stopping power in MeV * cm^2 / g
  const scalar dEdx{I.compute_bremsstrahlung(mat, ptc, {ptc, qop}) /
                    mat.mass_density() /
                    (unit<scalar>::MeV * unit<scalar>::cm2 / unit<scalar>::g)};

  const scalar expected_dEdx = std::get<3>(GetParam());

  // We have not implemented the bremsstrahlung for the heavier charged
  // particles, which is negligible
  if ((ptc.pdg_num() == electron<scalar>().pdg_num()) ||
      (ptc.pdg_num() == positron<scalar>().pdg_num())) {
    // Check if difference is within 14% error
    EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.14f);
  } else {
    EXPECT_FLOAT_EQ(static_cast<float>(dEdx), 0.f);
  }
}

// From https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
// Assumes that the stopping powers of electron and positron are the same

// Electrons
INSTANTIATE_TEST_SUITE_P(
    electron_Bremsstrahlung_100MeV_He, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.9229f)));

INSTANTIATE_TEST_SUITE_P(
    electron_Bremsstrahlung_100MeV_Al, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 3.714f)));

INSTANTIATE_TEST_SUITE_P(
    electron_Bremsstrahlung_100MeV_Si, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 4.099f)));

INSTANTIATE_TEST_SUITE_P(
    electron_Bremsstrahlung_100MeV_Cu, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), electron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 7.079f)));

// Positrons
INSTANTIATE_TEST_SUITE_P(
    positron_Bremsstrahlung_100MeV_He, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.9229f)));

INSTANTIATE_TEST_SUITE_P(
    positron_Bremsstrahlung_100MeV_Al, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 3.714f)));

INSTANTIATE_TEST_SUITE_P(
    positron_Bremsstrahlung_100MeV_Si, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 4.099f)));

INSTANTIATE_TEST_SUITE_P(
    positron_Bremsstrahlung_100MeV_Cu, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), positron<scalar>(),
                                      100.0f * unit<scalar>::MeV, 7.079f)));

// We have not implemented the bremsstrahlung for muons
INSTANTIATE_TEST_SUITE_P(
    muon_Bremsstrahlung_100MeV_Cu_muon, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.f)));

INSTANTIATE_TEST_SUITE_P(
    antimuon_Bremsstrahlung_100MeV_Cu_muon, EnergyLossBremsValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), antimuon<scalar>(),
                                      100.0f * unit<scalar>::MeV, 0.f)));
