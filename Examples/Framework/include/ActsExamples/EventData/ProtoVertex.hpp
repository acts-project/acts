// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Index.hpp"

#include <vector>

namespace ActsExamples {

/// A proto vertex is a collection of tracks identified by their indices.
using ProtoVertex = std::vector<Index>;
/// Container of proto vertices. Each proto vertex is identified by its index.
using ProtoVertexContainer = std::vector<ProtoVertex>;

}  // namespace ActsExamples
