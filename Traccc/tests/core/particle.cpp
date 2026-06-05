/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/utils/particle.hpp"

#include "traccc/definitions/common.hpp"
#include "traccc/edm/track_parameters.hpp"

// GTest include(s).
#include <gtest/gtest.h>

TEST(particle, particle_with_pdg_num) {
    const auto muon = traccc::detail::particle_from_pdg_number<float>(13);

    EXPECT_FLOAT_EQ(muon.mass(), 105.6583755f * traccc::unit<float>::MeV);
    EXPECT_FLOAT_EQ(muon.charge(), -1.f);

    const auto positron = traccc::detail::particle_from_pdg_number<float>(-11);

    EXPECT_FLOAT_EQ(positron.mass(), traccc::constant<float>::m_e);
    EXPECT_FLOAT_EQ(positron.charge(), 1.f);
}

TEST(particle, charge_conjugation) {

    const auto positron =
        traccc::detail::charge_conjugation(traccc::electron<float>());
    EXPECT_FLOAT_EQ(positron.mass(), traccc::constant<float>::m_e);
    EXPECT_FLOAT_EQ(positron.charge(), 1.f);
}

TEST(particle, correct_particle_hypothesis) {

    traccc::free_track_parameters<> free_trk({1.f, 1.f, 1.f}, 0.f,
                                             {2.f, 3.f, 4.f}, 1.f);

    const auto positron = traccc::detail::correct_particle_hypothesis(
        traccc::electron<float>(), free_trk);

    EXPECT_FLOAT_EQ(positron.mass(), traccc::constant<float>::m_e);
    EXPECT_FLOAT_EQ(positron.charge(), 1.f);
}
