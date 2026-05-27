/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/device/triplet_seeding_algorithm.hpp"

#include "traccc/edm/device/doublet_counter.hpp"
#include "traccc/edm/device/seeding_global_counter.hpp"

// Project include(s).
#include "traccc/seeding/spacepoint_binning_helper.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

struct triplet_seeding_algorithm::data {

    /// @name Configuration objects
    /// @{

    /// Configuration for the spacepoint grid forming step
    spacepoint_grid_config m_grid_config;
    /// Configuration for the seed finding step
    seedfinder_config m_finder_config;
    /// Configuration for the seed filtering step
    seedfilter_config m_filter_config;

    /// @}

    /// Axes for the internal spacepoint grid
    std::pair<details::spacepoint_grid_types::host::axis_p0_type,
              details::spacepoint_grid_types::host::axis_p1_type>
        m_axes;

};  // struct triplet_seeding_algorithm::data

triplet_seeding_algorithm::triplet_seeding_algorithm(
    const seedfinder_config& finder_config,
    const spacepoint_grid_config& grid_config,
    const seedfilter_config& filter_config, const memory_resource& mr,
    vecmem::copy& copy, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_data{std::make_unique<data>(
          grid_config, finder_config, filter_config,
          get_axes(grid_config, (mr.host ? *(mr.host) : mr.main)))} {}

triplet_seeding_algorithm::~triplet_seeding_algorithm() = default;

auto triplet_seeding_algorithm::operator()(
    const edm::spacepoint_collection::const_view& spacepoints) const
    -> output_type {

    // A small sanity check.
    assert(m_data);

    // Get the number of spacepoints. In an asynchronous way if possible.
    edm::spacepoint_collection::const_view::size_type n_spacepoints = 0u;
    if (mr().host) {
        vecmem::async_size size = copy().get_size(spacepoints, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_spacepoints = size.get();
    } else {
        n_spacepoints = copy().get_size(spacepoints);
    }

    // If there are no spacepoints, return right away.
    if (n_spacepoints == 0) {
        return {};
    }

    // Set up the container that will be filled with the required capacities for
    // the spacepoint grid.
    const unsigned int grid_bins =
        m_data->m_axes.first.n_bins * m_data->m_axes.second.n_bins;
    vecmem::data::vector_buffer<unsigned int> grid_capacities_buffer(grid_bins,
                                                                     mr().main);
    copy().setup(grid_capacities_buffer)->ignore();
    copy().memset(grid_capacities_buffer, 0)->ignore();

    // Launch the grid capacity counting kernel.
    count_grid_capacities_kernel({n_spacepoints, m_data->m_finder_config,
                                  m_data->m_axes.first, m_data->m_axes.second,
                                  spacepoints, grid_capacities_buffer});

    // Copy grid capacities back to the host
    vecmem::vector<unsigned int> grid_capacities_host(mr().host ? mr().host
                                                                : &(mr().main));
    copy()(grid_capacities_buffer, grid_capacities_host)->wait();

    // Create the spacepoint grid buffer and a prefix sum buffer that describes
    // it. (The latter is needed for the next few steps.)
    details::spacepoint_grid_types::buffer grid_buffer(
        m_data->m_axes.first, m_data->m_axes.second,
        std::vector<std::size_t>(grid_capacities_host.begin(),
                                 grid_capacities_host.end()),
        mr().main, mr().host, vecmem::data::buffer_type::resizable);
    copy().setup(grid_buffer._buffer)->ignore();
    vecmem::data::vector_buffer<prefix_sum_element_t> grid_prefix_sum_buffer(
        n_spacepoints, mr().main, vecmem::data::buffer_type::resizable);
    copy().setup(grid_prefix_sum_buffer)->ignore();

    // Launch the grid population kernel.
    populate_grid_kernel({n_spacepoints, m_data->m_finder_config, spacepoints,
                          grid_buffer, grid_prefix_sum_buffer});

    // Update the spacepoint counter. So that it would only count the "good"
    // spacepoints that ended up in the spacepoint grid.
    if (mr().host) {
        vecmem::async_size size =
            copy().get_size(grid_prefix_sum_buffer, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_spacepoints = size.get();
    } else {
        n_spacepoints = copy().get_size(grid_prefix_sum_buffer);
    }

    // Set up the doublet counter buffer.
    device::doublet_counter_collection_types::buffer doublet_counter_buffer{
        n_spacepoints, mr().main, vecmem::data::buffer_type::resizable};
    copy().setup(doublet_counter_buffer)->ignore();

    // Set up a global counter used in the seeding kernels.
    vecmem::unique_alloc_ptr<device::seeding_global_counter>
        globalCounter_device =
            vecmem::make_unique_alloc<device::seeding_global_counter>(
                mr().main);
    copy()
        .memset(vecmem::data::vector_view<device::seeding_global_counter>(
                    1u, globalCounter_device.get()),
                0)
        ->ignore();

    // Launch the doublet counting kernel.
    count_doublets_kernel(
        {n_spacepoints, m_data->m_finder_config, spacepoints, grid_buffer,
         grid_prefix_sum_buffer, doublet_counter_buffer,
         globalCounter_device->m_nMidBot, globalCounter_device->m_nMidTop});

    // Get the number of doublets found.
    device::doublet_counter_collection_types::buffer::size_type n_doublets = 0u;
    if (mr().host) {
        vecmem::async_size size =
            copy().get_size(doublet_counter_buffer, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_doublets = size.get();
    } else {
        n_doublets = copy().get_size(doublet_counter_buffer);
    }
    vecmem::unique_alloc_ptr<device::seeding_global_counter>
        globalCounter_host =
            vecmem::make_unique_alloc<device::seeding_global_counter>(
                mr().host ? *(mr().host) : mr().main);
    copy()(vecmem::data::vector_view<device::seeding_global_counter>(
               1u, globalCounter_device.get()),
           vecmem::data::vector_view<device::seeding_global_counter>(
               1u, globalCounter_host.get()))
        ->wait();

    // Exit already here if we won't find any triplets anyway.
    if ((globalCounter_host->m_nMidBot == 0) ||
        (globalCounter_host->m_nMidTop == 0)) {
        return {};
    }

    // Set up the doublet buffers.
    device_doublet_collection_types::buffer doublet_buffer_mb{
        globalCounter_host->m_nMidBot, mr().main};
    copy().setup(doublet_buffer_mb)->ignore();
    device_doublet_collection_types::buffer doublet_buffer_mt{
        globalCounter_host->m_nMidTop, mr().main};
    copy().setup(doublet_buffer_mt)->ignore();

    // Launch the doublet finding kernel.
    find_doublets_kernel({n_doublets, m_data->m_finder_config, spacepoints,
                          grid_buffer, doublet_counter_buffer,
                          doublet_buffer_mb, doublet_buffer_mt});

    // Set up the triplet counter buffers
    triplet_counter_spM_collection_types::buffer triplet_counter_spM_buffer{
        n_doublets, mr().main};
    copy().setup(triplet_counter_spM_buffer)->ignore();
    copy().memset(triplet_counter_spM_buffer, 0)->ignore();
    triplet_counter_collection_types::buffer triplet_counter_midBot_buffer{
        globalCounter_host->m_nMidBot, mr().main,
        vecmem::data::buffer_type::resizable};
    copy().setup(triplet_counter_midBot_buffer)->ignore();

    // Launch the triplet counting kernel.
    count_triplets_kernel({globalCounter_host->m_nMidBot,
                           m_data->m_finder_config, spacepoints, grid_buffer,
                           doublet_counter_buffer, doublet_buffer_mb,
                           doublet_buffer_mt, triplet_counter_spM_buffer,
                           triplet_counter_midBot_buffer});
    // Launch the triplet count reduction kernel.
    triplet_counts_reduction_kernel({n_doublets, doublet_counter_buffer,
                                     triplet_counter_spM_buffer,
                                     globalCounter_device->m_nTriplets});

    // Get the number of triplets found.
    triplet_counter_collection_types::buffer::size_type n_midBotTriplets = 0u;
    if (mr().host) {
        vecmem::async_size size =
            copy().get_size(triplet_counter_midBot_buffer, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        n_midBotTriplets = size.get();
    } else {
        n_midBotTriplets = copy().get_size(triplet_counter_midBot_buffer);
    }
    copy()(vecmem::data::vector_view<device::seeding_global_counter>(
               1u, globalCounter_device.get()),
           vecmem::data::vector_view<device::seeding_global_counter>(
               1u, globalCounter_host.get()))
        ->wait();

    // If no triplets could be found, exit already here.
    if (globalCounter_host->m_nTriplets == 0) {
        return {};
    }

    // Set up the triplet buffer.
    device_triplet_collection_types::buffer triplet_buffer{
        globalCounter_host->m_nTriplets, mr().main};
    copy().setup(triplet_buffer)->ignore();

    // Launch the triplet finding kernel.
    find_triplets_kernel({n_midBotTriplets, m_data->m_finder_config,
                          m_data->m_filter_config, spacepoints, grid_buffer,
                          doublet_counter_buffer, doublet_buffer_mt,
                          triplet_counter_spM_buffer,
                          triplet_counter_midBot_buffer, triplet_buffer});
    // Launch the triplet weight updating/filling kernel.
    update_triplet_weights_kernel(
        {globalCounter_host->m_nTriplets, m_data->m_filter_config, spacepoints,
         triplet_counter_spM_buffer, triplet_counter_midBot_buffer,
         triplet_buffer});

    // Create the result object.
    edm::seed_collection::buffer seed_buffer(
        globalCounter_host->m_nTriplets, mr().main,
        vecmem::data::buffer_type::resizable);
    copy().setup(seed_buffer)->ignore();

    // Launch the seed selecting/filling kernel.
    select_seeds_kernel(
        {n_doublets, m_data->m_finder_config, m_data->m_filter_config,
         spacepoints, grid_buffer, triplet_counter_spM_buffer,
         triplet_counter_midBot_buffer, triplet_buffer, seed_buffer});

    // Return the seed buffer.
    return seed_buffer;
}

}  // namespace traccc::device
