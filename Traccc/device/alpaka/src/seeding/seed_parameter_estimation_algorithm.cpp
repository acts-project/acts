/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
//
// Local include(s).
#include "traccc/alpaka/seeding/seed_parameter_estimation_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/seeding/device/estimate_track_params.hpp"

namespace traccc::alpaka {
namespace kernels {

/// Alpaka kernel for running @c traccc::device::estimate_track_params
template <typename bfield_t>
struct estimate_track_params {
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc, const track_params_estimation_config config,
        typename edm::measurement_collection::const_view measurements,
        edm::spacepoint_collection::const_view spacepoints,
        edm::seed_collection::const_view seeds, const bfield_t bfield,
        bound_track_parameters_collection_types::view params) const {
        auto const globalThreadIdx =
            ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Threads>(acc)[0u];

        device::estimate_track_params(globalThreadIdx, config, measurements,
                                      spacepoints, seeds, bfield, params);
    }
};  // struct estimate_track_params

}  // namespace kernels

seed_parameter_estimation_algorithm::seed_parameter_estimation_algorithm(
    const track_params_estimation_config& config,
    const traccc::memory_resource& mr, vecmem::copy& copy, alpaka::queue& q,
    std::unique_ptr<const Logger> logger)
    : device::seed_parameter_estimation_algorithm(config, mr, copy,
                                                  std::move(logger)),
      alpaka::algorithm_base(q) {}

void seed_parameter_estimation_algorithm::estimate_seed_params_kernel(
    const struct estimate_seed_params_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 4;
    const unsigned int n_blocks = (payload.n_seeds + n_threads - 1) / n_threads;
    magnetic_field_visitor<bfield_type_list<scalar>>(
        payload.bfield,
        [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
            ::alpaka::exec<Acc>(details::get_queue(queue()),
                                makeWorkDiv<Acc>(n_blocks, n_threads),
                                kernels::estimate_track_params<bfield_view_t>{},
                                payload.config, payload.measurements,
                                payload.spacepoints, payload.seeds, bfield,
                                payload.params);
        });
}

}  // namespace traccc::alpaka
