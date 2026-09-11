/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/csv/dfe.hpp"
#include "traccc/io/csv/hit.hpp"
#include "traccc/io/csv/measurement.hpp"
#include "traccc/io/csv/measurement_hit_id.hpp"
#include "traccc/io/csv/particle.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/simulation/measurement_smearer.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"

// Detray core include(s).
#include <detray/propagator/actors/parameter_updater.hpp>
#include <detray/propagator/base_actor.hpp>
#include <detray/utils/concepts.hpp>

// System include(s).
#include <filesystem>
#include <string>

namespace traccc {

template <typename smearer_t>
struct smearing_writer : detray::base_actor {
  using algebra_type = typename smearer_t::algebra_type;
  using scalar_type = detray::dscalar<algebra_type>;

  using measurement_hit_id_writer =
      io::dfe::NamedTupleCsvWriter<io::csv::measurement_hit_id>;
  using measurement_writer = io::dfe::NamedTupleCsvWriter<io::csv::measurement>;
  using hit_writer = io::dfe::NamedTupleCsvWriter<io::csv::hit>;
  using particle_writer = io::dfe::NamedTupleCsvWriter<io::csv::particle>;

  struct config {
    smearer_t smearer;
  };

  struct state {
    state(std::size_t event_id, const config& writer_cfg,
          const std::string directory)
        : m_particle_writer((std::filesystem::path{directory} /
                             traccc::io::get_event_filename(
                                 event_id, "-particles_initial.csv"))
                                .native()),
          m_hit_writer((std::filesystem::path{directory} /
                        traccc::io::get_event_filename(event_id, "-hits.csv"))
                           .native()),
          m_meas_writer(
              (std::filesystem::path{directory} /
               traccc::io::get_event_filename(event_id, "-measurements.csv"))
                  .native()),
          m_measurement_hit_id_writer(
              (std::filesystem::path{directory} /
               traccc::io::get_event_filename(event_id,
                                              "-measurement-simhit-map.csv"))
                  .native()),
          m_meas_smearer(writer_cfg.smearer) {}

    uint64_t particle_id = 0u;
    particle_writer m_particle_writer;
    hit_writer m_hit_writer;
    measurement_writer m_meas_writer;
    measurement_hit_id_writer m_measurement_hit_id_writer;
    uint64_t m_hit_count = 0u;
    smearer_t m_meas_smearer;

    void set_seed(const uint_fast64_t sd) { m_meas_smearer.set_seed(sd); }

    void write_particle(
        const traccc::free_track_parameters<algebra_type>& track,
        const traccc::pdg_particle<scalar_type>& ptc_type) {
      io::csv::particle particle;
      const auto pos = track.pos();
      const auto mom = track.mom(ptc_type.charge());

      particle.particle_id = particle_id;
      particle.particle_type = ptc_type.pdg_num();
      particle.vx = static_cast<float>(pos[0]);
      particle.vy = static_cast<float>(pos[1]);
      particle.vz = static_cast<float>(pos[2]);
      particle.vt = static_cast<float>(track.time());
      particle.px = static_cast<float>(mom[0]);
      particle.py = static_cast<float>(mom[1]);
      particle.pz = static_cast<float>(mom[2]);
      particle.q = static_cast<float>(ptc_type.charge());

      m_particle_writer.append(particle);
    }
  };

  struct measurement_kernel {
    template <typename mask_group_t, typename index_range_t>
    inline void operator()(
        const mask_group_t& /*mask_group*/, const index_range_t& /*index*/,
        const traccc::bound_track_parameters<algebra_type>& bound_params,
        smearer_t& smearer, io::csv::measurement& iomeas) const {
      using mask_t = typename mask_group_t::value_type;
      smearer.template operator()<mask_t>(smearer.get_offset(), bound_params,
                                          iomeas);
    }
  };

  template <typename propagator_state_t>
  void operator()(
      state& writer_state, propagator_state_t& propagation,
      const detray::actor::parameter_transporter_result<algebra_type>& res)
      const {
    auto& navigation = propagation.navigation();

    // triggered only for sensitive surfaces
    if (res.status <= detray::actor::status::e_success &&
        navigation.is_on_sensitive()) {
      const auto sf = navigation.current_surface();

      // Get the free param after the random scatterer update
      const auto& bound_params = res.destination_params();
      const auto free_param =
          sf.bound_to_free_vector(propagation.context(), bound_params);

      assert(sf.identifier() == bound_params.surface_link());

      TRACCC_DEBUG_HOST("Actor: Writing truth param.:\n" << bound_params);

      const auto pos = free_param.pos();
      const auto mom =
          free_param.mom(propagation.stepping().particle_hypothesis().charge());

      // Write hits
      io::csv::hit hit;
      hit.particle_id = writer_state.particle_id;
      hit.geometry_id = sf.identifier().value();
      hit.tx = static_cast<float>(pos[0]);
      hit.ty = static_cast<float>(pos[1]);
      hit.tz = static_cast<float>(pos[2]);
      hit.tt = static_cast<float>(free_param.time());
      hit.tpx = static_cast<float>(mom[0]);
      hit.tpy = static_cast<float>(mom[1]);
      hit.tpz = static_cast<float>(mom[2]);

      writer_state.m_hit_writer.append(hit);

      // Write measurements
      io::csv::measurement meas;

      assert(sf.identifier() == bound_params.surface_link());

      TRACCC_DEBUG_HOST("Actor: Writing truth param.:\n" << bound_params);

      meas.measurement_id = writer_state.m_hit_count;
      meas.geometry_id = hit.geometry_id;
      auto stddev_0 = writer_state.m_meas_smearer.stddev[0];
      auto stddev_1 = writer_state.m_meas_smearer.stddev[1];
      meas.var_local0 = static_cast<float>(stddev_0 * stddev_0);
      meas.var_local1 = static_cast<float>(stddev_1 * stddev_1);
      meas.phi = static_cast<float>(bound_params.phi());
      meas.theta = static_cast<float>(bound_params.theta());
      meas.time = static_cast<float>(bound_params.time());

      // Set local_key and smeared_local
      sf.template visit_mask<measurement_kernel>(
          bound_params, writer_state.m_meas_smearer, meas);

      writer_state.m_meas_writer.append(meas);

      // Write hit measurement map
      io::csv::measurement_hit_id measurement_hit_id;
      measurement_hit_id.hit_id = writer_state.m_hit_count;
      measurement_hit_id.measurement_id = writer_state.m_hit_count;
      writer_state.m_measurement_hit_id_writer.append(measurement_hit_id);
      writer_state.m_hit_count++;
    }
  }
};
}  // namespace traccc
