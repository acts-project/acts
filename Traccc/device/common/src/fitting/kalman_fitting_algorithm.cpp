/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/fitting/device/kalman_fitting_algorithm.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_buffer.hpp>

// System include(s).
#include <numeric>
#include <vector>

namespace traccc::device {

struct kalman_fitting_algorithm::data {

    /// @name Configuration object(s)
    /// @{

    /// Configuration for the fitting algorithm
    fitting_config m_config;

    /// @}

};  // struct kalman_fitting_algorithm::data

kalman_fitting_algorithm::fit_payload::fit_payload(const detector_buffer& det,
                                                   const magnetic_field& f)
    : detector(det), field(f) {}

kalman_fitting_algorithm::kalman_fitting_algorithm(
    const config_type& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_data{std::make_unique<data>(config)} {}

kalman_fitting_algorithm::~kalman_fitting_algorithm() = default;

kalman_fitting_algorithm::output_type kalman_fitting_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& field,
    const edm::track_container<default_algebra>::const_view& input_tracks)
    const {

    // Get the number of tracks and the number of constituens (measurements)
    // in each track. Note that the number of tracks is not "resizable". That
    // we can just get directly from the view object. But the number of
    // constituents need to be taken from the buffer itself.
    assert(copy().get_size(input_tracks.tracks) ==
           input_tracks.tracks.capacity());
    const edm::track_collection<default_algebra>::const_view::size_type
        n_tracks = input_tracks.tracks.capacity();
    if (n_tracks == 0) {
        // Return early, if there are no tracks.
        return {};
    }
    std::vector<unsigned int> candidate_sizes;
    if (mr().host) {
        vecmem::async_sizes sizes =
            copy().get_sizes(input_tracks.tracks, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        auto& temp = sizes.get();
        candidate_sizes = {temp.begin(), temp.end()};
    } else {
        candidate_sizes = copy().get_sizes(input_tracks.tracks);
    }

    // Get the total number of states (measurements) to fit.
    const unsigned int n_states =
        std::accumulate(candidate_sizes.begin(), candidate_sizes.end(), 0u);

    // Create the result buffer.
    edm::track_container<default_algebra>::buffer output_tracks{
        {candidate_sizes, mr().main, mr().host,
         vecmem::data::buffer_type::resizable},
        {n_states, mr().main, vecmem::data::buffer_type::resizable},
        input_tracks.measurements};
    copy().setup(output_tracks.tracks)->ignore();
    copy().setup(output_tracks.states)->ignore();

    // Create the order to fit the tracks in.
    vecmem::data::vector_buffer<device::sort_key> track_sort_keys(n_tracks,
                                                                  mr().main);
    vecmem::data::vector_buffer<unsigned int> track_indices{n_tracks,
                                                            mr().main};
    prepare_track_fit_order(input_tracks.tracks, track_sort_keys,
                            track_indices);

    // Create the buffer(s) used during the fitting.
    vecmem::data::vector_buffer<unsigned int> track_liveness(n_tracks,
                                                             mr().main);

    // Run "fitting prelude" kernel.
    fit_prelude_kernel(
        {track_indices, input_tracks, output_tracks, track_liveness});

    // Calculate the number of surfaces to use during the fit for each track.
    // Then create the concrete buffer such that it could be passed to the
    // fitting functions through a polymorphic pointer.
    std::vector<unsigned int> n_surfaces(candidate_sizes.size());
    std::transform(candidate_sizes.begin(), candidate_sizes.end(),
                   n_surfaces.begin(), [&](const unsigned int sz) {
                       return std::max(
                           sz * m_data->m_config.surface_sequence_size_factor,
                           m_data->m_config.min_surface_sequence_capacity);
                   });

    // Prepare the payload for the fitting kernel(s).
    const fit_payload payload = prepare_fit_payload(
        det, field, n_surfaces, {track_indices, track_liveness, output_tracks});

    // Run the iterative track fitting.
    for (std::size_t i = 0; i < m_data->m_config.n_iterations; ++i) {
        fit_forward_kernel(m_data->m_config, payload);
        fit_backward_kernel(m_data->m_config, payload);
    }

    // Return the fitted tracks.
    return output_tracks;
}

}  // namespace traccc::device
