/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/simulation/smearing_writer.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"
#include "traccc/utils/propagation.hpp"

// Detray include(s).
#include <detray/test/utils/random_scatterer.hpp>

// System include(s).
#include <limits>
#include <memory>

namespace traccc {

template <typename detector_t, typename bfield_t, typename track_generator_t,
          typename writer_t>
struct simulator {

    using algebra_type = typename detector_t::algebra_type;
    using scalar_type = typename detector_t::scalar_type;
    using bfield_type = bfield_t;

    struct config {
        detray::propagation::config propagation;

        /// Particle hypothesis
        traccc::pdg_particle<scalar_type> ptc_type{traccc::muon<scalar_type>()};

        // Simulation setup
        bool do_energy_loss = true;
        bool do_multiple_scattering = true;
        bool m_is_min_pT = false;
        scalar_type m_min_p = 10.f * traccc::unit<scalar_type>::MeV;

        /// Set the momentum limit to @param p
        inline void min_p(const scalar_type p) {
            m_is_min_pT = false;
            m_min_p = p;
        }

        /// Set the transverse momentum limit to @param p
        inline void min_pT(const scalar_type p) {
            m_is_min_pT = true;
            m_min_p = p;
        }
    };

    using actor_chain_type = detray::actor_chain<
        detray::actor::momentum_aborter<scalar_type>,
        detray::actor::parameter_updater<
            algebra_type, detray::actor::random_scatterer<algebra_type>,
            writer_t>>;

    using navigator_type = detray::caching_navigator<detector_t>;
    using stepper_type = detray::rk_stepper<
        typename bfield_type::view_t, algebra_type,
        detray::constrained_step<detray::dscalar<algebra_type>>>;
    using propagator_type =
        detray::propagator<stepper_type, navigator_type, actor_chain_type>;

    simulator(const detray::pdg_particle<scalar>& ptc_type, std::size_t events,
              const detector_t& det, const bfield_type& field,
              track_generator_t&& track_gen,
              typename writer_t::config&& writer_cfg,
              const std::string directory = "")
        : m_events(events),
          m_directory(directory),
          m_detector(det),
          m_field(field),
          m_track_generator(
              std::make_unique<track_generator_t>(std::move(track_gen))),
          m_writer_cfg(writer_cfg) {

        m_cfg.ptc_type = ptc_type;
        m_track_generator->config().charge(ptc_type.charge());

        // Turn off tracking features
        m_cfg.propagation.stepping.do_covariance_transport = false;
        m_cfg.propagation.stepping.use_eloss_gradient = false;
        m_cfg.propagation.stepping.use_field_gradient = false;
        m_cfg.propagation.navigation.estimate_scattering_noise = false;

        m_updater_state = detray::actor::parameter_updater_state<algebra_type>{
            m_cfg.propagation};
    }

    config& get_config() { return m_cfg; }

    void run() {

        TRACCC_VERBOSE_HOST("Running fast simulation...");

        // Update the actor config
        if (m_cfg.m_is_min_pT) {
            m_aborter_state.min_pT(m_cfg.m_min_p);
            m_aborter_state.min_p(0.f);
        } else {
            m_aborter_state.min_p(m_cfg.m_min_p);
            m_aborter_state.min_pT(0.f);
        }
        m_scatterer.do_energy_loss = m_cfg.do_energy_loss;
        m_scatterer.do_multiple_scattering = m_cfg.do_multiple_scattering;

        for (std::size_t event_id = 0u; event_id < m_events; event_id++) {

            typename writer_t::state writer_state(
                event_id, std::move(m_writer_cfg), m_directory);

            // Set random seed
            m_scatterer.set_seed(event_id);
            writer_state.set_seed(event_id);

            auto actor_states = detray::tie(m_aborter_state, m_updater_state,
                                            m_scatterer, writer_state);

            for (auto track : *m_track_generator.get()) {

                m_updater_state.init(track);

                writer_state.write_particle(
                    track,
                    detail::correct_particle_hypothesis(m_cfg.ptc_type, track));

                typename propagator_type::state propagation(track, m_field,
                                                            m_detector);
                propagation.set_particle(
                    detail::correct_particle_hypothesis(m_cfg.ptc_type, track));

                propagator_type p(m_cfg.propagation);

                // Set overstep tolerance and stepper constraint
                propagation.stepping()
                    .template set_constraint<
                        detray::step::constraint::e_accuracy>(
                        m_cfg.propagation.stepping.step_constraint);

                p.propagate(propagation, actor_states);

                // Increase the particle id
                writer_state.particle_id++;
            }
        }
    }

    private:
    config m_cfg;
    std::size_t m_events{0u};
    std::string m_directory = "";
    const detector_t& m_detector;
    const typename bfield_type::view_t m_field;
    std::unique_ptr<track_generator_t> m_track_generator;
    typename writer_t::config m_writer_cfg;

    /// Actor states
    typename detray::actor::momentum_aborter<scalar_type>::state
        m_aborter_state{};
    typename detray::actor::parameter_updater<algebra_type>::state
        m_updater_state{};
    typename detray::actor::random_scatterer<algebra_type>::state m_scatterer{};
};

}  // namespace traccc
