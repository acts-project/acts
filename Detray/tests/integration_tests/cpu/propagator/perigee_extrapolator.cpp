/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/propagator/perigee_extrapolator.hpp"

#include "detray/definitions/units.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

using namespace detray;

using test_algebra = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;

constexpr scalar tol{1e-3f};

namespace detray::test {

struct portal_aborter : public base_actor {
    /// Exits the navigation as soon as the first portal is reached.
    /// @param prop_state state of the propagation
    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(propagator_state_t &prop_state) const {
        DETRAY_VERBOSE_HOST("Aborter: Check target surface");

        auto &navigation = prop_state.navigation();
        if (navigation.is_on_portal()) {
            navigation.exit();
            prop_state.heartbeat(false);
        }
    }
};

}  // namespace detray::test

/// Test the extrapolation of a track to the perigee
GTEST_TEST(detray_propagator, perigee_extrapolator) {

    vecmem::host_memory_resource host_mr;

    toy_det_config<scalar> toy_cfg{};
    toy_cfg.use_material_maps(false);
    const auto [d, names] = build_toy_detector<test_algebra>(host_mr, toy_cfg);

    using detector_t = decltype(d);
    using stepper_t = line_stepper<test_algebra>;
    using actor_chain_t = actor_chain<test::portal_aborter>;

    using propagator_t =
        propagator<stepper_t, caching_navigator<detector_t>, actor_chain_t>;
    using extrapolator_t = perigee_extrapolator<detector_t, stepper_t>;

    // Track generator configuration
    using generator_t =
        uniform_track_generator<free_track_parameters<test_algebra>>;

    // Track generator config
    generator_t::configuration trk_gen_cfg{};
    trk_gen_cfg.eta_range(-3.f, 3.f);
    trk_gen_cfg.phi_steps(50u).eta_steps(50u);
    trk_gen_cfg.p_tot(10.f * unit<scalar>::GeV);

    propagation::config prop_cfg{};

    propagator_t p{prop_cfg};
    extrapolator_t pe{prop_cfg};

    // Iterate through uniformly distributed momentum directions
    for (auto track : generator_t{trk_gen_cfg}) {
        // Propagate to first portal
        propagator_t::state prop_state{track, d, prop_cfg.context};
        ASSERT_TRUE(p.propagate(prop_state));
        ASSERT_TRUE(prop_state.stepping().path_length() >=
                    25.f * unit<scalar>::mm);
        ASSERT_NEAR(vector::perp(prop_state.stepping()().pos()),
                25.f * unit<scalar>::mm, tol);

        // Propagate back to perigee
        extrapolator_t::state extrapolator_state =
            pe.create_state(prop_state.stepping()(), d);

        // Was the minimum distance to the first portal covered?
        ASSERT_NEAR(vector::perp(extrapolator_state.stepping()().pos()),
                25.f * unit<scalar>::mm, tol);

        // Extrapolation success?
        const auto bound_params = pe.extrapolate(extrapolator_state);
        ASSERT_TRUE(pe.finished(extrapolator_state));
        // Was the minimum distance to the first portal covered?
        EXPECT_TRUE(extrapolator_state.stepping().abs_path_length() >=
                    25.f * unit<scalar>::mm);
        // Back at origin?
        EXPECT_NEAR(vector::norm(extrapolator_state.stepping()().pos() -
                                 point3{0.f, 0.f, 0.f}),
                0.f, tol);

        const scalar d0{bound_params[e_bound_loc0]};
        const scalar z0{bound_params[e_bound_loc1]};

        EXPECT_NEAR(d0, 0.f, tol);
        EXPECT_NEAR(z0, 0.f, tol);
    }
}
