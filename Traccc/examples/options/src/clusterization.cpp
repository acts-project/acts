/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/clusterization.hpp"

#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/examples/utils/printable.hpp"

namespace traccc::opts {

clusterization::clusterization() : interface("Clusterization Options") {

    m_desc.add_options()(
        "threads-per-partition",
        boost::program_options::value(&m_config.threads_per_partition)
            ->default_value(m_config.threads_per_partition),
        "The number of threads per partition");
    m_desc.add_options()(
        "max-cells-per-thread",
        boost::program_options::value(&m_config.max_cells_per_thread)
            ->default_value(m_config.max_cells_per_thread),
        "The maximum number of cells per thread");
    m_desc.add_options()(
        "target-cells-per-thread",
        boost::program_options::value(&m_config.target_cells_per_thread)
            ->default_value(m_config.target_cells_per_thread),
        "The target number of cells per thread");
    m_desc.add_options()(
        "backup-size-multiplier",
        boost::program_options::value(&m_config.backup_size_multiplier)
            ->default_value(m_config.backup_size_multiplier),
        "The size multiplier of the backup scratch space");
}

clusterization::operator clustering_config() const {
    return m_config;
}

clusterization::operator host::clusterization_algorithm::config_type() const {
    return {};
}

std::unique_ptr<configuration_printable> clusterization::as_printable() const {
    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Threads per partition",
        std::to_string(m_config.threads_per_partition)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Target cells per thread",
        std::to_string(m_config.target_cells_per_thread)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Max cells per thread", std::to_string(m_config.max_cells_per_thread)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Scratch space multiplier",
        std::to_string(m_config.backup_size_multiplier)));

    return cat;
}
}  // namespace traccc::opts
