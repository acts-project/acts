// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
