/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/clusterization/device/clusterization_algorithm.hpp"

namespace traccc::alpaka {

/// Algorithm performing hit clusterization
///
/// This algorithm implements hit clusterization in a massively-parallel
/// approach. Each thread handles a pre-determined number of detector cells.
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying this buffer.
///
class clusterization_algorithm : public device::clusterization_algorithm,
                                 public algorithm_base {

    public:
    /// Constructor for clusterization algorithm
    ///
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param q The Alpaka queue to perform the operations in
    /// @param config The clustering configuration
    ///
    clusterization_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy, alpaka::queue& q,
        const config_type& config,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    private:
    /// @name Function(s) inherited from the base class
    /// @{

    /// Function meant to perform sanity checks on the input data
    ///
    /// @param cells     All cells in an event
    /// @return @c true if the input data is valid, @c false otherwise
    ///
    bool input_is_valid(
        const edm::silicon_cell_collection::const_view& cells) const override;

    /// Main CCL kernel launcher
    void ccl_kernel(const ccl_kernel_payload& payload) const override;

    /// Cluster data reification kernel launcher
    ///
    /// @param num_cells    Number of cells in the event
    /// @param disjoint_set Buffer for the disjoint set data structure
    /// @param cluster_data The cluster collection to fill
    ///
    void cluster_maker_kernel(
        unsigned int num_cells,
        const vecmem::data::vector_view<unsigned int>& disjoint_set,
        edm::silicon_cluster_collection::view& cluster_data) const override;

    /// @}

};  // class clusterization_algorithm

}  // namespace traccc::alpaka
