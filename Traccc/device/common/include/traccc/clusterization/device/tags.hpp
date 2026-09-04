/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {
/**
 * @defgroup Clustering algorithm cluster retention parameters
 *
 * Optional parameters to clustering algorithms which convert directly from
 * cells to measurements, determining whether to reconstruct the intermediate
 * cluster data or not.
 *
 * @{
 * @brief Explicitly discard cluster information.
 */
struct clustering_discard_disjoint_set {};
/*
 * @brief Explicitly reconstruct and return cluster information.
 */
struct clustering_keep_disjoint_set {};
/*
 * @}
 */
}  // namespace traccc::device
