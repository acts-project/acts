// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

namespace ActsPython {
namespace Concepts {
template <typename T>
concept has_write_method =
    requires(T& t, const ActsExamples::AlgorithmContext& ctx) {
      { t.write(ctx) } -> std::same_as<ActsExamples::ProcessCode>;
    };
}  // namespace Concepts
}  // namespace ActsPython