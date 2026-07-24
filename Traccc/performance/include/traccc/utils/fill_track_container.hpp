/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/data_format.hpp"
#include "traccc/utils/event_data.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/transcribe_to_trace.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <vector>

namespace traccc {

/// Fill track containers Kalman comparison from simulation data
void fill_track_containers(
    std::unique_ptr<const traccc::Logger> ilogger,
    const traccc::host_detector* host_det, const std::string& input_dir,
    const unsigned int n_events, const bool use_acts_geoid,
    const traccc::scalar track_min_pT, const traccc::scalar track_max_rad,
    std::vector<traccc::free_track_parameters<traccc::default_algebra>>& tracks,
    std::vector<vecmem::vector<traccc::propagation_validator::candidate_type<
        traccc::default_detector::host>>>& truth_traces,
    typename edm::measurement_collection::host& measurements,
    traccc::edm::track_container<traccc::default_algebra>::host&
        track_container) {
  using algebra_t = traccc::default_algebra;
  using detector_t = traccc::default_detector::host;

  TRACCC_LOCAL_LOGGER(std::move(ilogger));

  // Host memory resource
  vecmem::host_memory_resource host_mr;

  // Geometry context
  const traccc::default_detector::host::geometry_context ctx{};

  // Retrieve detector
  const detector_t& det = host_det->template as<traccc::default_detector>();

  tracks.reserve(n_events * 1000);
  truth_traces.reserve(tracks.capacity());

  TRACCC_VERBOSE("Reconstructing " << n_events << " events");

  // Read the truth data in
  for (std::size_t i_event = 0u; i_event < n_events; ++i_event) {
    traccc::event_data evt_data(input_dir, i_event, host_mr, use_acts_geoid,
                                host_det, data_format::csv);
    TRACCC_VERBOSE("Event " << i_event << ": Found "
                            << evt_data.m_particle_map.size()
                            << " initial particles");

    if (evt_data.m_particle_map.empty()) {
      TRACCC_ERROR("Removing event " << i_event << ": Found no particles");
      continue;
    }

    if (evt_data.m_ptc_to_meas_map.empty()) {
      TRACCC_ERROR(
          "Removing event "
          << i_event
          << ": Found no connections between particles and measurements");
      continue;
    }

    std::size_t n_tracks{tracks.size()};
    for (const auto& [ptc_id, ptc] : evt_data.m_particle_map) {
      // Particle produced no valid measurements
      if (!evt_data.m_ptc_to_meas_map.contains(ptc)) {
        TRACCC_WARNING("Event "
                       << i_event << ": Removing particle " << ptc_id
                       << " since it has no measurements linked to it");
        continue;
      }
      // Minimum momentum
      const traccc::scalar pT{vector::perp(ptc.momentum)};
      if (pT <= track_min_pT) {
        TRACCC_WARNING("Event " << i_event << ": Removing particle " << ptc_id
                                << " due to transv. momentum cut (pT was "
                                << pT / traccc::unit<traccc::scalar>::MeV
                                << " MeV)");
        continue;
      }

      // Make a trace of detray-understandable intersections
      auto truth_trace_fw = traccc::propagation_validator::transcribe_to_trace(
          ctx, det, ptc, evt_data.m_ptc_to_meas_map);

      // No meansurements found in event
      if (truth_trace_fw.empty()) {
        TRACCC_WARNING("Event " << i_event
                                << ": No measurements found for particle "
                                << ptc_id);
        continue;
      }

      // Minimum radius (remove secondaries)
      const traccc::scalar rad{vector::perp(ptc.vertex)};
      if (rad >= track_max_rad) {
        TRACCC_WARNING("Event " << i_event << ": Removing particle " << ptc_id
                                << " due to radius cut (radius was "
                                << rad / traccc::unit<traccc::scalar>::mm
                                << " mm)");
        continue;
      }

      assert(!truth_trace_fw.empty());

      truth_traces.push_back(std::move(truth_trace_fw));

      TRACCC_DEBUG("Event " << i_event << ": Found " << truth_traces.size()
                            << " truth measurement(s) for track "
                            << tracks.size());

      // Transcribe measurements to global collection
      const auto& measurements_per_ptc = evt_data.m_ptc_to_meas_map.at(ptc);
      assert(!measurements_per_ptc.empty());

      // Current min measurement index for this track
      const auto meas_offset{static_cast<unsigned int>(measurements.size())};

      for (const auto& meas : measurements_per_ptc) {
        measurements.push_back(meas);
      }

      // Device measurement container with the correct(!) size
      edm::measurement_collection::const_device device_measurements{
          vecmem::get_data(measurements)};

      // Create the track and link one track state per measurement to it
      edm::track_collection<algebra_t>::host::object_type track;
      for (const auto& [i, meas] :
           detray::views::enumerate(measurements_per_ptc)) {
        const auto meas_idx{static_cast<unsigned int>(i + meas_offset)};
        const auto state_idx{
            static_cast<unsigned int>(track_container.states.size())};

        track.constituent_links().push_back(
            {edm::track_constituent_link::track_state, state_idx});

        track_container.states.push_back(
            edm::make_track_state<algebra_t>(device_measurements, meas_idx));
      }

      // Record the truth track candidate.
      track_container.tracks.push_back(track);

      TRACCC_DEBUG("-> Truth tracks " << track_container.tracks.size());

      // Construct the initial free track parameters for the propagation
      tracks.emplace_back(ptc.vertex, 0.f, ptc.momentum, ptc.charge);
    }

    if (tracks.size() - n_tracks == 0u) {
      TRACCC_WARNING("Event " << i_event << ": No eligible tracks in event");
    } else {
      TRACCC_VERBOSE("Event " << i_event << ": Found "
                              << tracks.size() - n_tracks
                              << " reconstructible truth track(s) in event");
    }
  }

  // Make the final measurement info available to the input containers
  track_container.measurements = vecmem::get_data(measurements);
}

}  // namespace traccc
