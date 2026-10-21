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
#include "detray/material/material_slab.hpp"
#include "detray/material/predefined_materials.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/landau_distribution.hpp"
#include "detray/test/utils/random_scatterer.hpp"
#include "detray/test/utils/statistics.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <random>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;

// Test class for energy loss with Bethe function
// Input tuple: < material, particle type, momentum, expected output >
// https://pdg.lbl.gov/2022/AtomicNuclearProperties for muon dEdX and range
class EnergyLossBetheValidation
    : public ::testing::TestWithParam<
          std::tuple<material<scalar>, pdg_particle<scalar>, scalar, scalar>> {
};

// This tests the material functionalities
TEST_P(EnergyLossBetheValidation, bethe_energy_loss) {
  // Interaction object
  interaction<scalar> I;

  // Zero incidence angle
  const scalar cos_inc_ang{1.f};

  // Material slab with the 1 cm thickness
  material_slab<scalar> slab(std::get<0>(GetParam()), 1.f * unit<scalar>::cm);

  // Particle
  pdg_particle<scalar> ptc = std::get<1>(GetParam());

  // Path segment in the material
  const scalar path_segment{slab.path_segment(cos_inc_ang)};

  // qOverP
  const scalar qop{ptc.charge() / std::get<2>(GetParam())};

  // Bethe Stopping power in MeV * cm^2 / g
  const scalar dEdx{I.compute_energy_loss_bethe_bloch(
                        path_segment, slab.get_material(), ptc, {ptc, qop}) /
                    path_segment / slab.get_material().mass_density() /
                    (unit<scalar>::MeV * unit<scalar>::cm2 / unit<scalar>::g)};

  const scalar expected_dEdx = std::get<3>(GetParam());

  // Check if difference is within 5% error
  EXPECT_NEAR((expected_dEdx - dEdx) / dEdx, 0.f, 0.05f);
}

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 6.539f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 4.182f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 4.777f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_H2Liquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(hydrogen_liquid<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 5.305f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 3.082f)));

/*
//@fixme: Test fails with He Gas and 1 GeV muons (18 % difference)
INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 2.133f)));
*/

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.768f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_HeGas, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(helium_gas<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 3.188f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.533f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.744f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.097f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Al, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(aluminium<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.360f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.608f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.803f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.177f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Si, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.451f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_0p1GeV_Si_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             silicon_with_ded<scalar>(), muon<scalar>(),
                             0.1003f * unit<scalar>::GeV, 2.608f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_1GeV_Si_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             silicon_with_ded<scalar>(), muon<scalar>(),
                             1.101f * unit<scalar>::GeV, 1.803f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_10GeV_Si_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             silicon_with_ded<scalar>(), muon<scalar>(),
                             10.11f * unit<scalar>::GeV, 2.177f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_100GeV_Si_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             silicon_with_ded<scalar>(), muon<scalar>(),
                             100.1f * unit<scalar>::GeV, 2.451f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_ArLiquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(argon_liquid<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.34f)));

// @fixme ~6% discrepancy
/*INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_ArLiquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(argon_liquid<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.644f)));
*/

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_ArLiquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(argon_liquid<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 2.003f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_ArLiquid, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(argon_liquid<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.258f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Fe, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.274f)));

// @fixme ~6% discrepancy
/*INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Fe, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.581f)));
*/

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Fe, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 1.942f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Fe, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.207f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Fe_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron_with_ded<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.274f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Fe_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron_with_ded<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.581f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Fe_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron_with_ded<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 1.942f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Fe_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(iron_with_ded<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.207f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Cu, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.198f)));

// @fixme
/*INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Cu, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.532f)));
*/

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Cu, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 1.891f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Cu, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.155f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_Cu_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper_with_ded<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 2.198f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_Cu_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper_with_ded<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.532f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_Cu_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper_with_ded<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 1.891f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_Cu_with_DED, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(copper_with_ded<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.155f)));

// @fixme ~10% discrepancy
/*INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_0p1GeV_CsI, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(cesium_iodide<scalar>(), muon<scalar>(),
                                      0.1003f * unit<scalar>::GeV, 1.869f)));
*/

// @fixme ~10% discrepancy
/*INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_1GeV_CsI, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(cesium_iodide<scalar>(), muon<scalar>(),
                                      1.101f * unit<scalar>::GeV, 1.391f)));
*/

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_10GeV_CsI, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(cesium_iodide<scalar>(), muon<scalar>(),
                                      10.11f * unit<scalar>::GeV, 1.755f)));

INSTANTIATE_TEST_SUITE_P(
    detray_material_Bethe_100GeV_CsI, EnergyLossBetheValidation,
    ::testing::Values(std::make_tuple(cesium_iodide<scalar>(), muon<scalar>(),
                                      100.1f * unit<scalar>::GeV, 2.012f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_0p1GeV_CsI_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             cesium_iodide_with_ded<scalar>(), muon<scalar>(),
                             0.1003f * unit<scalar>::GeV, 1.869f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_1GeV_CsI_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             cesium_iodide_with_ded<scalar>(), muon<scalar>(),
                             1.101f * unit<scalar>::GeV, 1.391f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_10GeV_CsI_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             cesium_iodide_with_ded<scalar>(), muon<scalar>(),
                             10.11f * unit<scalar>::GeV, 1.755f)));

INSTANTIATE_TEST_SUITE_P(detray_material_Bethe_100GeV_CsI_with_DED,
                         EnergyLossBetheValidation,
                         ::testing::Values(std::make_tuple(
                             cesium_iodide_with_ded<scalar>(), muon<scalar>(),
                             100.1f * unit<scalar>::GeV, 2.012f)));

// Test class for MUON energy loss with Landau function
// Input tuple: < material / particle type / energy / expected energy loss  /
// expected fwhm  >
class EnergyLossLandauValidation
    : public ::testing::TestWithParam<std::tuple<
          material<scalar>, pdg_particle<scalar>, scalar, scalar, scalar>> {};

TEST_P(EnergyLossLandauValidation, landau_energy_loss) {
  // Interaction object
  interaction<scalar> I;

  // Zero incidence angle
  const scalar cos_inc_ang{1.f};

  // Material
  const auto mat = std::get<0>(GetParam());

  // Thickness
  const scalar thickness = 0.17f * unit<scalar>::cm;

  // Material slab with a unit thickness
  material_slab<scalar> slab(mat, thickness);

  // Particle
  pdg_particle<scalar> ptc = std::get<1>(GetParam());

  // Path segment in the material
  const scalar path_segment{slab.path_segment(cos_inc_ang)};

  // Energy
  const scalar E = std::get<2>(GetParam());

  // p
  const scalar p = std::sqrt(E * E - ptc.mass() * ptc.mass());

  // q
  const scalar q = ptc.charge();

  // qOverP
  const scalar qop{q / p};

  // Landau Energy loss in MeV
  const scalar Landau_MeV{I.compute_energy_loss_landau(path_segment,
                                                       slab.get_material(), ptc,
                                                       {ptc, qop}) /
                          unit<scalar>::MeV};

  // Check if difference is within 5% error
  EXPECT_TRUE(std::abs(std::get<3>(GetParam()) - Landau_MeV) / Landau_MeV <
              0.05f);

  // Landau Energy loss Fluctuation
  const scalar fwhm_MeV{I.compute_energy_loss_landau_fwhm(path_segment,
                                                          slab.get_material(),
                                                          ptc, {ptc, qop}) /
                        unit<scalar>::MeV};

  // Check if difference is within 10% error
  EXPECT_TRUE(std::abs(std::get<4>(GetParam()) - fwhm_MeV) / fwhm_MeV < 0.1f);
}

// Expected output from Fig 33.7 in RPP2018
INSTANTIATE_TEST_SUITE_P(
    detray_material_Landau_10GeV_Silicon, EnergyLossLandauValidation,
    ::testing::Values(std::make_tuple(silicon<scalar>(), muon<scalar>(),
                                      10.f * unit<scalar>::GeV, 0.525f,
                                      0.13f)));

// Input tuple: < energy >
class LandauDistributionValidation
    : public ::testing::TestWithParam<
          std::tuple<pdg_particle<scalar>, scalar>> {};

TEST_P(LandauDistributionValidation, landau_distribution) {
  // Random generator
  std::random_device rd{};
  std::mt19937_64 generator{rd()};
  generator.seed(0u);

  // Interaction object
  interaction<scalar> I;

  // Zero incidence angle
  const scalar cos_inc_ang{1.f};

  // Material
  const auto mat = silicon<scalar>();

  // Thickness
  const auto thickness = 1.f * unit<scalar>::mm;

  // Material slab with a unit thickness
  material_slab<scalar> slab(mat, thickness);

  // Particle
  pdg_particle<scalar> ptc = std::get<0>(GetParam());

  // Path segment in the material
  const scalar path_segment{slab.path_segment(cos_inc_ang)};

  // Energy
  const scalar E = std::get<1>(GetParam());

  // p
  const scalar p = std::sqrt(E * E - ptc.mass() * ptc.mass());

  // q
  const scalar q = ptc.charge();

  // qOverP
  const scalar qop{q / p};

  // Bethe energy loss
  const scalar dE{I.compute_energy_loss_bethe_bloch(
      path_segment, slab.get_material(), ptc, {ptc, qop})};

  // Landau Energy loss
  const scalar mpv{I.compute_energy_loss_landau(
      path_segment, slab.get_material(), ptc, {ptc, qop})};

  // Landau Energy loss Sigma
  const scalar sigma{I.compute_energy_loss_landau_sigma(
      path_segment, slab.get_material(), ptc, {ptc, qop})};

  // Landau Sampling
  landau_distribution<scalar> landau;

  std::vector<scalar> samples;
  std::size_t n_samples{10000u};
  for (std::size_t i = 0u; i < n_samples; i++) {
    const scalar sa = landau(generator, mpv, sigma);
    samples.push_back(sa);
  }

  // Mean energy loss from Landau samples
  const scalar sample_mean = statistics::mean(samples);

  // Make sure that the difference is within 50% (...)
  // @note: We should not expect a good consistency between landau
  // distribution samples and bethe bloch function because the arguments
  // for landau_distribution are not the most probable value and sigma
  // (sigma is not even defined in Landau distribution). We might need to put
  // a smart scaling factor
  EXPECT_NEAR((dE - sample_mean) / dE, 0.f, 0.5f);

  // Test attenuate function
  std::vector<scalar> energies;

  const scalar mass = ptc.mass();

  for (std::size_t i = 0u; i < n_samples; i++) {
    const scalar new_p = actor::random_scatterer<test_algebra>().attenuate(
        mpv, sigma, mass, p, generator);
    ASSERT_TRUE(new_p < p);

    const scalar new_E = std::sqrt(new_p * new_p + mass * mass);
    energies.push_back(new_E);
  }

  const scalar E_mean = statistics::mean(energies);
  const scalar E_expected = E - dE;
  EXPECT_TRUE(E_expected < E);

  // Make sure that the difference is within 65%
  EXPECT_NEAR((E_mean - E_expected) / dE, 0.f, 0.65f);
}

INSTANTIATE_TEST_SUITE_P(
    detray_material_Landau, LandauDistributionValidation,
    ::testing::Values(std::make_tuple(muon<scalar>(), 1.f * unit<scalar>::GeV),
                      std::make_tuple(muon<scalar>(), 10.f * unit<scalar>::GeV),
                      std::make_tuple(muon<scalar>(),
                                      100.f * unit<scalar>::GeV)));
