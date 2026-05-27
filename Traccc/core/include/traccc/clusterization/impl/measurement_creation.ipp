/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/utils/detray_conversion.hpp"

// System include(s).
#include <cassert>

namespace traccc::details {

template <typename TCell, typename TDesign>
TRACCC_HOST_DEVICE inline vector2 position_from_cell(
    const edm::silicon_cell<TCell>& cell,
    const traccc::detector_design_description_interface<TDesign>& module_dd) {
    // Calculate / construct the local cell position.
    vector2 cell_lower_position = {
        (module_dd.bin_edges_x()).at(cell.channel0()),
        (module_dd.bin_edges_y()).at(cell.channel1())};

    vector2 cell_upper_position = {
        (module_dd.bin_edges_x()).at(cell.channel0() + 1),
        (module_dd.bin_edges_y()).at(cell.channel1() + 1)};

    vector2 cell_middle_position = {
        scalar{0.5f} * (cell_upper_position[0] + cell_lower_position[0]),
        scalar{0.5f} * (cell_upper_position[1] + cell_lower_position[1])};

    return cell_middle_position;
}

template <typename T, typename TDesign>
TRACCC_HOST_DEVICE inline void calc_cluster_properties(
    const edm::silicon_cluster<T>& cluster,
    const edm::silicon_cell_collection::const_device& cells,
    const traccc::detector_design_description_interface<TDesign>& module_dd,
    point2& mean, point2& var, scalar& totalWeight) {

    point2 offset{0.f, 0.f}, width{0.f, 0.f};
    bool first_processed = false;

    unsigned int min_channel0 = std::numeric_limits<unsigned int>::max();
    unsigned int max_channel0 = std::numeric_limits<unsigned int>::lowest();
    unsigned int min_channel1 = std::numeric_limits<unsigned int>::max();
    unsigned int max_channel1 = std::numeric_limits<unsigned int>::lowest();

    // Loop over the cell indices of the cluster.
    for (const unsigned int cell_idx : cluster.cell_indices()) {

        // The cell object.
        const edm::silicon_cell cell = cells.at(cell_idx);

        // Translate the cell readout value into a weight.
        const scalar weight = cell.activation();

        // Update all output properties with this cell.
        totalWeight += weight;
        scalar weight_factor = weight / totalWeight;

        point2 cell_position = position_from_cell(cell, module_dd);

        min_channel0 = std::min(min_channel0, cell.channel0());
        min_channel1 = std::min(min_channel1, cell.channel1());
        max_channel0 = std::max(max_channel0, cell.channel0());
        max_channel1 = std::max(max_channel1, cell.channel1());

        if (!first_processed) {
            offset = cell_position;
            first_processed = true;
        }

        cell_position = cell_position - offset;

        const point2 diff_old = cell_position - mean;
        mean = mean + diff_old * weight_factor;
        const point2 diff_new = cell_position - mean;

        var[0] = (1.f - weight_factor) * var[0] +
                 weight_factor * (diff_old[0] * diff_new[0]);
        var[1] = (1.f - weight_factor) * var[1] +
                 weight_factor * (diff_old[1] * diff_new[1]);
    }

    // cluster width in the number of cells
    unsigned int delta0 = (max_channel0 - min_channel0) + 1;
    unsigned int delta1 = (max_channel1 - min_channel1) + 1;

    vector2 cluster_lower_position = {
        (module_dd.bin_edges_x()).at(min_channel0),
        (module_dd.bin_edges_y()).at(min_channel1)};

    vector2 cluster_upper_position = {
        (module_dd.bin_edges_x()).at(max_channel0 + 1),
        (module_dd.bin_edges_y()).at(max_channel1 + 1)};

    width[0] = cluster_upper_position[0] - cluster_lower_position[0];
    width[1] = cluster_upper_position[1] - cluster_lower_position[1];

    point2 pitch = {width[0] / static_cast<float>(delta0),
                    width[1] / static_cast<float>(delta1)};

    var = var + point2{pitch[0] * pitch[0] / static_cast<scalar>(12.),
                       pitch[1] * pitch[1] / static_cast<scalar>(12.)};

    mean = mean + offset;
}

template <typename T1, typename T2>
TRACCC_HOST_DEVICE inline void fill_measurement(
    edm::measurement<T1>& measurement, const edm::silicon_cluster<T2>& cluster,
    const unsigned int index,
    const edm::silicon_cell_collection::const_device& cells,
    const detector_design_description::const_device& det_descr,
    const detector_conditions_description::const_device& det_cond) {

    // To calculate the mean and variance with high numerical stability
    // we use a weighted variant of Welford's algorithm. This is a
    // single-pass online algorithm that works well for large numbers
    // of samples, as well as samples with very high values.
    //
    // To learn more about this algorithm please refer to:
    // [1] https://doi.org/10.1080/00401706.1962.10490022
    // [2] The Art of Computer Programming, Donald E. Knuth, second
    //     edition, chapter 4.2.2.

    // Security checks.
    assert(cluster.cell_indices().empty() == false);
    assert([&]() {
        const unsigned int module_idx =
            cells.module_index().at(cluster.cell_indices().front());
        for (const unsigned int cell_idx : cluster.cell_indices()) {
            if (cells.module_index().at(cell_idx) != module_idx) {
                return false;
            }
        }
        return true;
    }() == true);

    // The index of the module the cluster is on.
    const unsigned int module_idx =
        cells.module_index().at(cluster.cell_indices().front());
    const auto module_cd = det_cond.at(module_idx);
    const unsigned int design_idx = module_cd.module_to_design_id();
    const auto module_dd = det_descr.at(design_idx);

    // Calculate the cluster properties
    scalar totalWeight = 0.f;
    point2 mean{0.f, 0.f}, var{0.f, 0.f};
    calc_cluster_properties(cluster, cells, module_dd, mean, var, totalWeight);
    assert(totalWeight > 0.f);

    // Fill the measurement object.
    measurement.surface_link() = module_cd.geometry_id();

    // apply lorentz shift to the cell position
    measurement.local_position() = utils::to_float_array<default_algebra>(
        mean + module_cd.measurement_translation());

    // plus pitch^2 / 12
    measurement.local_variance() = utils::to_float_array<default_algebra>(var);

    // For the ambiguity resolution algorithm, give a unique measurement ID
    measurement.identifier() = index;
    measurement.cluster_index() = index;

    // Set the measurement dimensionality.
    measurement.dimensions() = module_dd.dimensions();

    // Set the measurement's subspace.
    measurement.set_subspace(module_dd.subspace());

    // Save the index of the cluster that produced this measurement
    measurement.cluster_index() = static_cast<unsigned int>(index);
}

}  // namespace traccc::details
