/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/clusterization/device/ccl_kernel_definitions.hpp"
#include "traccc/clusterization/device/tags.hpp"
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/unique_ptr.hpp>
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>
#include <optional>

namespace traccc::device {

/// Base class for the algorithms performing hit clusterization
///
/// This algorithm implements hit clusterization in a massively-parallel
/// approach. Each thread handles a pre-determined number of detector cells.
///
/// This algorithm returns a buffer which is not necessarily filled yet. A
/// synchronisation statement is required before destroying the buffer.
///
class clusterization_algorithm
    : public algorithm<edm::measurement_collection::buffer(
          const edm::silicon_cell_collection::const_view&,
          const detector_design_description::const_view&,
          const detector_conditions_description::const_view&)>,
      public algorithm<edm::measurement_collection::buffer(
          const edm::silicon_cell_collection::const_view&,
          const detector_design_description::const_view&,
          const detector_conditions_description::const_view&,
          clustering_discard_disjoint_set&&)>,
      public algorithm<std::pair<edm::measurement_collection::buffer,
                                 edm::silicon_cluster_collection::buffer>(
          const edm::silicon_cell_collection::const_view&,
          const detector_design_description::const_view&,
          const detector_conditions_description::const_view&,
          clustering_keep_disjoint_set&&)>,
      public messaging,
      public algorithm_base {

    public:
    /// Configuration type
    using config_type = clustering_config;

    /// Constructor for clusterization algorithm
    ///
    /// @param mr The memory resource(s) to use in the algorithm
    /// @param copy The copy object to use for copying data between device
    ///             and host memory blocks
    /// @param config The clustering configuration
    /// partition
    ///
    clusterization_algorithm(
        const traccc::memory_resource& mr, vecmem::copy& copy,
        const config_type& config,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Callable operator for clusterization algorithm
    ///
    /// @param cells     All cells in an event
    /// @param det_descr The detector segmnentation description
    /// @param det_cond The detector conditions description
    /// @return a measurement collection (buffer)
    ///
    /// @{
    edm::measurement_collection::buffer operator()(
        const edm::silicon_cell_collection::const_view& cells,
        const detector_design_description::const_view& det_descr,
        const detector_conditions_description::const_view& det_cond)
        const override;

    edm::measurement_collection::buffer operator()(
        const edm::silicon_cell_collection::const_view& cells,
        const detector_design_description::const_view& det_descr,
        const detector_conditions_description::const_view& det_cond,
        clustering_discard_disjoint_set&&) const override;

    std::pair<edm::measurement_collection::buffer,
              edm::silicon_cluster_collection::buffer>
    operator()(const edm::silicon_cell_collection::const_view& cells,
               const detector_design_description::const_view& det_descr,
               const detector_conditions_description::const_view& det_cond,
               clustering_keep_disjoint_set&&) const override;
    /// @}

    protected:
    /// @name Function(s) to be implemented by derived classes
    /// @{

    /// Function meant to perform sanity checks on the input data
    ///
    /// @param cells     All cells in an event
    /// @return @c true if the input data is valid, @c false otherwise
    ///
    virtual bool input_is_valid(
        const edm::silicon_cell_collection::const_view& cells) const = 0;

    /// Payload for the @c ccl_kernel function
    struct ccl_kernel_payload {
        /// Number of cells in the event
        unsigned int n_cells;
        /// The clustering configuration
        const config_type& config;
        /// All cells in an event
        const edm::silicon_cell_collection::const_view& cells;
        /// The detector description
        const detector_design_description::const_view& det_descr;
        /// The detector conditions description
        const detector_conditions_description::const_view& det_cond;
        /// The measurement collection to fill
        edm::measurement_collection::view& measurements;
        /// Buffer for backup of the first element links
        vecmem::data::vector_view<details::fallback_index_t>& f_backup;
        /// Buffer for backup of the group first element links
        vecmem::data::vector_view<details::fallback_index_t>& gf_backup;
        /// Buffer for backup of the adjacency matrix (counts)
        vecmem::data::vector_view<unsigned char>& adjc_backup;
        /// Buffer for backup of the adjacency matrix (values)
        vecmem::data::vector_view<details::fallback_index_t>& adjv_backup;
        /// Mutex for the backup structures
        unsigned int* backup_mutex;
        /// Buffer for the disjoint set data structure
        vecmem::data::vector_view<unsigned int>& disjoint_set;
        /// Buffer for the sizes of the clusters
        vecmem::data::vector_view<unsigned int>& cluster_sizes;
    };

    /// Main CCL kernel launcher
    ///
    /// @param payload The payload containing all necessary data for the kernel
    ///
    virtual void ccl_kernel(const ccl_kernel_payload& payload) const = 0;

    /// Cluster data reification kernel launcher
    ///
    /// @param num_cells    Number of cells in the event
    /// @param disjoint_set Buffer for the disjoint set data structure
    /// @param cluster_data The cluster collection to fill
    ///
    virtual void cluster_maker_kernel(
        unsigned int num_cells,
        const vecmem::data::vector_view<unsigned int>& disjoint_set,
        edm::silicon_cluster_collection::view& cluster_data) const = 0;

    /// @}

    private:
    /// Main algorithmic implementation of the clusterization algorithm
    std::pair<edm::measurement_collection::buffer,
              std::optional<edm::silicon_cluster_collection::buffer>>
    execute_impl(const edm::silicon_cell_collection::const_view& cells,
                 const detector_design_description::const_view& det_descr,
                 const detector_conditions_description::const_view& det_cond,
                 bool keep_disjoint_set) const;

    /// Clusterization configuration
    config_type m_config;
    /// Memory reserved for edge cases
    mutable vecmem::data::vector_buffer<details::fallback_index_t> m_f_backup;
    mutable vecmem::data::vector_buffer<details::fallback_index_t> m_gf_backup;
    mutable vecmem::unique_alloc_ptr<unsigned int> m_backup_mutex;
    mutable vecmem::data::vector_buffer<unsigned char> m_adjc_backup;
    mutable vecmem::data::vector_buffer<details::fallback_index_t>
        m_adjv_backup;

};  // class clusterization_algorithm

}  // namespace traccc::device
