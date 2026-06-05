/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/resolution/res_plot_tool_config.hpp"
#include "traccc/resolution/stat_plot_tool_config.hpp"
#include "traccc/utils/messaging.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/particle.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_collection.hpp"
#include "traccc/utils/event_data.hpp"

// System include(s).
#include <cassert>
#include <memory>

namespace traccc {
namespace details {

/// Data members that should not pollute the API of
/// @c traccc::fitting_performance_writer
struct fitting_performance_writer_data;

}  // namespace details

class fitting_performance_writer : public messaging {

    public:
    struct config {
        /// Output filename.
        std::string file_path = "performance_track_fitting.root";
        /// Output file mode
        std::string file_mode = "RECREATE";
        /// Plot tool configurations.
        res_plot_tool_config res_config;
        stat_plot_tool_config stat_config;
    };

    /// Constructor with writer config
    fitting_performance_writer(const config& cfg,
                               std::unique_ptr<const traccc::Logger> logger);

    /// Destructor that closes the file
    ~fitting_performance_writer();

    /// Fill the tracking results into the histograms
    ///
    /// @param track The fitted track to measure/write the performance of
    /// @param track_states All reconstructed track states
    /// @param measurements All reconstructed measurements
    /// @param det detector object
    /// @param evt_map event map to find the truth values
    template <typename detector_t>
    void write(const edm::track_collection<
                   traccc::default_algebra>::host::proxy_type track,
               const edm::track_state_collection<traccc::default_algebra>::host&
                   track_states,
               const edm::measurement_collection::host& measurements,
               const detector_t& det, event_data& evt_data,
               const detector_t::geometry_context& ctx = {}) {

        static_assert(std::same_as<typename detector_t::algebra_type,
                                   traccc::default_algebra>);

        if (track.fit_outcome() != track_fit_outcome::SUCCESS) {
            return;
        }

        // Get the first smoothed track state
        const unsigned int trk_state_idx =
            std::find_if(track.constituent_links().begin(),
                         track.constituent_links().end(),
                         [&](const edm::track_constituent_link& link) {
                             assert(link.type ==
                                    edm::track_constituent_link::track_state);
                             return track_states.at(link.index).is_smoothed();
                         })
                ->index;
        const edm::track_state trk_state = track_states.at(trk_state_idx);
        assert(!trk_state.is_hole());
        assert(trk_state.is_smoothed());

        bool use_found = !evt_data.m_found_meas_to_ptc_map.empty();

        const std::map<event_data::measurement_proxy,
                       std::map<particle, std::size_t>>& meas_to_ptc_map =
            use_found ? evt_data.m_found_meas_to_ptc_map
                      : evt_data.m_meas_to_ptc_map;
        const std::map<event_data::measurement_proxy,
                       std::pair<point3, point3>>& meas_to_param_map =
            use_found ? evt_data.m_found_meas_to_param_map
                      : evt_data.m_meas_to_param_map;

        const edm::measurement meas =
            measurements.at(trk_state.measurement_index());

        // Find the contributing particle
        // @todo: Use identify_contributing_particles function
        std::map<particle, std::size_t> contributing_particles =
            meas_to_ptc_map.at(meas);

        const particle ptc = contributing_particles.begin()->first;

        // Find the truth global position and momentum
        const point3 global_pos = meas_to_param_map.at(meas).first;
        const vector3 global_mom = meas_to_param_map.at(meas).second;

        const detray::tracking_surface sf{det, meas.surface_link()};
        const point2 truth_bound =
            sf.global_to_bound(ctx, global_pos, vector::normalize(global_mom));

        // Return value
        bound_track_parameters<> truth_param{};
        truth_param.set_bound_local(truth_bound);
        truth_param.set_phi(vector::phi(global_mom));
        truth_param.set_theta(vector::theta(global_mom));
        // @todo: Assign a proper value to time
        truth_param.set_time(0.f);
        truth_param.set_qop(ptc.charge / vector::norm(global_mom));

        // For the moment, only fill with the first measurements
        write_res(truth_param, trk_state.smoothed_params(), ptc);
        write_stat(track, track_states, measurements);
    }

    /// Writing caches into the file
    void finalize();

    private:
    /// Non-templated part of the @c write(...) function
    void write_res(const bound_track_parameters<>& truth_param,
                   const bound_track_parameters<>& fit_param,
                   const particle& ptc);

    /// Non-templated part of the @c write(...) function
    void write_stat(
        const edm::track_collection<traccc::default_algebra>::host::proxy_type
            track,
        const edm::track_state_collection<traccc::default_algebra>::host&
            track_states,
        const edm::measurement_collection::host& measurements);

    /// Configuration for the tool
    config m_cfg;

    /// Opaque data members for the class
    std::unique_ptr<details::fitting_performance_writer_data> m_data;

};  // class fitting_performance_writer

}  // namespace traccc
