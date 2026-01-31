// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Hashing/AnnoyForwardDeclarations.hpp"

#include <cstdint>

namespace ActsPlugins {
/// @addtogroup hashing_plugin
/// @{

using AnnoyMetric = Annoy::AngularEuclidean;
using AnnoyModel =
    Annoy::AnnoyIndex<std::uint32_t, float, AnnoyMetric, Annoy::Kiss32Random,
                      Annoy::AnnoyIndexSingleThreadedBuildPolicy>;

/// @}
}  // namespace ActsPlugins
