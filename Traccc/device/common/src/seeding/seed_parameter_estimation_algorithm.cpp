/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/device/seed_parameter_estimation_algorithm.hpp"

namespace traccc::device {

struct seed_parameter_estimation_algorithm::data {

    /// Configuration for the track parameter estimation step
    track_params_estimation_config m_config;

};  // struct seed_parameter_estimation_algorithm::data

seed_parameter_estimation_algorithm::seed_parameter_estimation_algorithm(
    const track_params_estimation_config& config,
    const traccc::memory_resource& mr, vecmem::copy& copy,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_data{std::make_unique<data>(data{config})} {}

seed_parameter_estimation_algorithm::~seed_parameter_estimation_algorithm() =
    default;

auto seed_parameter_estimation_algorithm::operator()(
    const magnetic_field& bfield,
    const edm::measurement_collection::const_view& measurements,
    const edm::spacepoint_collection::const_view& spacepoints,
    const edm::seed_collection::const_view& seeds) const -> output_type {

    // Get the number of seeds. In an asynchronous way if possible.
    edm::seed_collection::const_view::size_type n_seeds = 0u;
    if (mr().host) {
        const vecmem::async_size size = copy().get_size(seeds, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_seeds = size.get();
    } else {
        n_seeds = copy().get_size(seeds);
    }

    // If there are no seeds, return right away.
    if (n_seeds == 0) {
        return {};
    }

    // Set up the output buffer.
    bound_track_parameters_collection_types::buffer result(n_seeds, mr().main);
    copy().setup(result)->ignore();

    // Launch the seed parameter estimation kernel.
    estimate_seed_params_kernel({n_seeds, m_data->m_config, bfield,
                                 measurements, spacepoints, seeds, result});

    // Return the result.
    return result;
}

}  // namespace traccc::device
