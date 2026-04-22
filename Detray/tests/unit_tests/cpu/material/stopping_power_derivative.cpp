// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
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

using scalar = test::scalar;

GTEST_TEST(detray_material, derivative_test_beta2) {
  const pdg_particle<scalar> ptc = muon<scalar>();

  // mass
  const scalar m = ptc.mass();

  // charge
  const scalar q = ptc.charge();

  // Displacement for numerical differentiaion
  const scalar h = 0.001f;

  // Iterate from 1 GeV to 10 GeV
  for (unsigned int i = 1u; i < 10; i++) {
    const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
    const scalar qop = q / p;

    detail::relativistic_quantities rq(m, qop, q);
    detail::relativistic_quantities rq1(m, qop + h, q);
    detail::relativistic_quantities rq2(m, qop - h, q);

    const scalar numerical = (rq1.m_beta2 - rq2.m_beta2) / (2.f * h);

    const scalar evaluated = rq.derive_beta2();

    EXPECT_NEAR((numerical - evaluated) / numerical, 0, 0.012);
  }
}

// Test class for the derivative of bethe stopping power
class DerivativeOfBetheEquationValidation
    : public ::testing::TestWithParam<
          std::tuple<pdg_particle<scalar>, material<scalar>>> {};

TEST_P(DerivativeOfBetheEquationValidation,
       derivative_of_bethe_stopping_power) {
  // Interaction object
  interaction<scalar> Interactor;

  const pdg_particle<scalar> ptc = std::get<0>(GetParam());

  // Material
  const auto mat = std::get<1>(GetParam());

  // mass
  const scalar m = ptc.mass();

  // charge
  const scalar q = ptc.charge();

  // Displacement for numerical differentiaion
  const scalar h = 1e-3f;

  // Iterate from 2 GeV to 100 GeV
  for (unsigned int i = 2u; i < 100; i++) {
    const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
    const scalar qop = q / p;

    const detail::relativistic_quantities<scalar> rq(m, qop, q);
    const detail::relativistic_quantities<scalar> rq1(m, qop + h, q);
    const detail::relativistic_quantities<scalar> rq2(m, qop - h, q);

    // Log term
    const scalar log_term1 = rq1.compute_bethe_bloch_log_term(mat);
    const scalar log_term2 = rq2.compute_bethe_bloch_log_term(mat);

    const scalar numerical_log_term = (log_term1 - log_term2) / (2.f * h);
    const scalar evaluated_log_term = rq.derive_bethe_bloch_log_term();

    EXPECT_NEAR(numerical_log_term, evaluated_log_term,
                numerical_log_term * 0.01f);

    // delta half
    const scalar dhalf1 = rq1.compute_delta_half(mat);
    const scalar dhalf2 = rq2.compute_delta_half(mat);

    const scalar numerical_dhalf = (dhalf1 - dhalf2) / (2.f * h);
    const scalar evaluated_dhalf = rq.derive_delta_half(mat);

    EXPECT_NEAR(numerical_dhalf, evaluated_dhalf, numerical_dhalf * 0.01f);

    // Bethe equation
    const scalar bethe1 = Interactor.compute_bethe_bloch(mat, ptc, rq1);
    const scalar bethe2 = Interactor.compute_bethe_bloch(mat, ptc, rq2);

    const scalar numerical_bethe = (bethe1 - bethe2) / (2.f * h);

    const scalar evaluated_bethe = Interactor.derive_bethe_bloch(mat, ptc, rq);

    EXPECT_NEAR(numerical_bethe, evaluated_bethe, numerical_bethe * 0.01f);
  }
}

INSTANTIATE_TEST_SUITE_P(
    detray_material_BetheEquationDerivative,
    DerivativeOfBetheEquationValidation,
    ::testing::Values(
        std::make_tuple(muon<scalar>(), hydrogen_gas<scalar>()),
        std::make_tuple(muon<scalar>(), helium_gas<scalar>()),
        std::make_tuple(muon<scalar>(), isobutane<scalar>()),
        std::make_tuple(muon<scalar>(), aluminium<scalar>()),
        std::make_tuple(muon<scalar>(), silicon<scalar>()),
        std::make_tuple(muon<scalar>(), tungsten<scalar>()),
        std::make_tuple(muon<scalar>(), gold<scalar>()),
        std::make_tuple(muon<scalar>(), cesium_iodide<scalar>()),
        std::make_tuple(muon<scalar>(), silicon_with_ded<scalar>()),
        std::make_tuple(muon<scalar>(), cesium_iodide_with_ded<scalar>())));

// Test class for the derivative of bremsstrahlung stopping power
class DerivativeOfBremsstrahlungValidation
    : public ::testing::TestWithParam<
          std::tuple<pdg_particle<scalar>, material<scalar>>> {};

TEST_P(DerivativeOfBremsstrahlungValidation,
       derivative_of_brems_stopping_power) {
  // Interaction object
  interaction<scalar> Interactor;

  const pdg_particle<scalar> ptc = std::get<0>(GetParam());

  // Material
  const auto mat = std::get<1>(GetParam());

  // mass
  const scalar m = ptc.mass();

  // charge
  const scalar q = ptc.charge();

  // Displacement for numerical differentiaion
  const scalar h = 1e-3f;

  // Iterate from 1 GeV to 100 GeV
  for (unsigned int i = 1u; i < 100; i++) {
    const scalar p = static_cast<scalar>(i) * detray::unit<scalar>::GeV;
    const scalar qop = q / p;

    const detail::relativistic_quantities<scalar> rq(m, qop, q);
    const detail::relativistic_quantities<scalar> rq1(m, qop + h, q);
    const detail::relativistic_quantities<scalar> rq2(m, qop - h, q);

    // Brems equation
    const scalar brems1 = Interactor.compute_bremsstrahlung(mat, ptc, rq1);
    const scalar brems2 = Interactor.compute_bremsstrahlung(mat, ptc, rq2);

    const scalar numerical_brems = (brems1 - brems2) / (2.f * h);

    const scalar evaluated_brems =
        Interactor.derive_bremsstrahlung(mat, ptc, rq);

    EXPECT_NEAR(numerical_brems, evaluated_brems, numerical_brems * 0.01f);
  }
}

INSTANTIATE_TEST_SUITE_P(
    detray_material_BremsstrahlungDerivative,
    DerivativeOfBremsstrahlungValidation,
    ::testing::Values(
        std::make_tuple(electron<scalar>(), hydrogen_gas<scalar>()),
        std::make_tuple(electron<scalar>(), helium_gas<scalar>()),
        std::make_tuple(electron<scalar>(), isobutane<scalar>()),
        std::make_tuple(electron<scalar>(), aluminium<scalar>()),
        std::make_tuple(electron<scalar>(), silicon<scalar>()),
        std::make_tuple(electron<scalar>(), tungsten<scalar>()),
        std::make_tuple(electron<scalar>(), gold<scalar>()),
        std::make_tuple(electron<scalar>(), cesium_iodide<scalar>()),
        std::make_tuple(electron<scalar>(), silicon_with_ded<scalar>()),
        std::make_tuple(electron<scalar>(), cesium_iodide_with_ded<scalar>())));
