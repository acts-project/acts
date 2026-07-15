/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/clusterization/device/clusterization_algorithm.hpp"

#include <limits>
#include <stdexcept>

namespace traccc::device {

clusterization_algorithm::clusterization_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& cp,
    const config_type& config, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, cp},
      m_config{config},
      m_f_backup{m_config.backup_size(), mr.main},
      m_gf_backup{m_config.backup_size(), mr.main},
      m_backup_mutex{vecmem::make_unique_alloc<unsigned int>(mr.main)},
      m_adjc_backup{m_config.backup_size(), mr.main},
      m_adjv_backup{m_config.backup_size() * 4, mr.main} {

    if (config.max_partition_size() >
        static_cast<unsigned int>(
            std::numeric_limits<details::index_t>::max())) {
        throw std::domain_error(
            "Maximal partition size is too large for the chosen index type!");
    }

    copy().setup(m_f_backup)->wait();
    copy().setup(m_gf_backup)->wait();
    copy().setup(m_adjc_backup)->wait();
    copy().setup(m_adjv_backup)->wait();
    copy()
        .memset(
            vecmem::data::vector_view<unsigned int>{1, m_backup_mutex.get()}, 0)
        ->wait();
}

edm::measurement_collection::buffer clusterization_algorithm::operator()(
    const edm::silicon_cell_collection::const_view& cells,
    const detector_design_description::const_view& det_descr,
    const detector_conditions_description::const_view& det_cond) const {

    return this->operator()(cells, det_descr, det_cond,
                            clustering_discard_disjoint_set{});
}

edm::measurement_collection::buffer clusterization_algorithm::operator()(
    const edm::silicon_cell_collection::const_view& cells,
    const detector_design_description::const_view& det_descr,
    const detector_conditions_description::const_view& det_cond,
    clustering_discard_disjoint_set&&) const {

    static constexpr bool KEEP_DISJOINT_SET = false;
    auto [res, djs] =
        this->execute_impl(cells, det_descr, det_cond, KEEP_DISJOINT_SET);
    assert(!djs.has_value());
    return std::move(res);
}

std::pair<edm::measurement_collection::buffer,
          edm::silicon_cluster_collection::buffer>
clusterization_algorithm::operator()(
    const edm::silicon_cell_collection::const_view& cells,
    const detector_design_description::const_view& det_descr,
    const detector_conditions_description::const_view& det_cond,
    clustering_keep_disjoint_set&&) const {

    static constexpr bool KEEP_DISJOINT_SET = true;
    auto [res, djs] =
        this->execute_impl(cells, det_descr, det_cond, KEEP_DISJOINT_SET);
    assert(djs.has_value());
    return {std::move(res), std::move(*djs)};
}

std::pair<edm::measurement_collection::buffer,
          std::optional<edm::silicon_cluster_collection::buffer>>
clusterization_algorithm::execute_impl(
    const edm::silicon_cell_collection::const_view& cells,
    const detector_design_description::const_view& det_descr,
    const detector_conditions_description::const_view& det_cond,
    bool keep_disjoint_set) const {

    // Check the input data in debug mode.
    assert(input_is_contiguous(cells));

    // Get the number of cells, in an asynchronous way if possible.
    edm::silicon_cell_collection::const_view::size_type num_cells = 0u;
    if (mr().host) {
        const vecmem::async_size size = copy().get_size(cells, *(mr().host));
        // Here we could give control back to the caller, once our code allows
        // for it. (coroutines...)
        num_cells = size.get();
    } else {
        num_cells = copy().get_size(cells);
    }

    // If there are no cells, return right away.
    if (num_cells == 0) {
        if (keep_disjoint_set) {
            return {edm::measurement_collection::buffer{},
                    edm::silicon_cluster_collection::buffer{}};
        } else {
            return {};
        }
    }

    // Create the result object, overestimating the number of measurements.
    edm::measurement_collection::buffer measurements{
        num_cells, mr().main, vecmem::data::buffer_type::resizable};
    copy().setup(measurements)->ignore();

    // Ensure that the chosen maximum cell count is compatible with the maximum
    // stack size.
    assert(m_config.max_cells_per_thread <=
           device::details::CELLS_PER_THREAD_STACK_LIMIT);

    // If we are keeping the disjoint set data structure, allocate space for it.
    vecmem::data::vector_buffer<unsigned int> disjoint_set;
    vecmem::data::vector_buffer<unsigned int> cluster_sizes;
    if (keep_disjoint_set) {
        disjoint_set = {num_cells, mr().main};
        cluster_sizes = {num_cells, mr().main};
    }

    std::optional<edm::silicon_cell_collection::buffer> sorted_cells;
    vecmem::data::vector_buffer<unsigned int> permutation_map_buffer;
    edm::silicon_cell_collection::const_view sorted_cells_view;

    if (m_config.sort_cells) {
        sorted_cells =
            edm::silicon_cell_collection::buffer(num_cells, mr().main);
        copy().setup(*sorted_cells)->ignore();
        // The permutation map contents are only needed to reify the cluster
        // data, but the buffer is allocated unconditionally because
        // implementations may use it as scratch space while sorting.
        permutation_map_buffer =
            vecmem::data::vector_buffer<unsigned int>(num_cells, mr().main);
        copy().setup(permutation_map_buffer)->ignore();

        sort_cells_kernel(num_cells, cells, *sorted_cells,
                          permutation_map_buffer, keep_disjoint_set);
        sorted_cells_view =
            edm::silicon_cell_collection::const_view(*sorted_cells);
    } else {
        sorted_cells_view = edm::silicon_cell_collection::const_view(cells);
    }

    // Sorting must preserve module contiguity; a violation here indicates a
    // defect in the sorting kernel (e.g. key construction), which the
    // ordering check alone cannot detect.
    assert(input_is_contiguous(sorted_cells_view));
    assert(input_is_sorted(sorted_cells_view));

    // Launch the CCL kernel.
    ccl_kernel({num_cells, m_config, sorted_cells_view, det_descr, det_cond,
                measurements, m_f_backup, m_gf_backup, m_adjc_backup,
                m_adjv_backup, m_backup_mutex.get(), disjoint_set,
                cluster_sizes});

    std::optional<edm::silicon_cluster_collection::buffer> cluster_data =
        std::nullopt;

    // Create the cluster data if requested.
    if (keep_disjoint_set) {

        // Get the number of reconstructed measurements, in an asynchronous way
        // if possible.
        edm::measurement_collection::buffer::size_type num_measurements = 0u;
        if (mr().host) {
            const vecmem::async_size size =
                copy().get_size(measurements, *(mr().host));
            // Here we could give control back to the caller, once our code
            // allows for it. (coroutines...)
            num_measurements = size.get();
        } else {
            num_measurements = copy().get_size(measurements);
        }

        // This could be further optimized by only copying the number of
        // elements necessary. But since cluster making is mainly meant for
        // performance measurements, on first order this should be good enough.
        vecmem::vector<unsigned int> cluster_sizes_host =
            ((mr().host != nullptr) ? vecmem::vector<unsigned int>(mr().host)
                                    : vecmem::vector<unsigned int>());
        copy()(cluster_sizes, cluster_sizes_host)->wait();
        cluster_sizes_host.resize(num_measurements);

        // Create the result cluster collection.
        cluster_data.emplace(cluster_sizes_host, mr().main, mr().host,
                             vecmem::data::buffer_type::resizable);
        copy().setup(*cluster_data)->ignore();

        // Run the cluster data reification kernel.
        cluster_maker_kernel(num_cells, disjoint_set, *cluster_data,
                             permutation_map_buffer);
    }

    // Return the reconstructed measurements.
    return {std::move(measurements), std::move(cluster_data)};
}

}  // namespace traccc::device
