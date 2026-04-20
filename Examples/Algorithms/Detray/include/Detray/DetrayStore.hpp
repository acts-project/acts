// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <memory>

#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace ActsExamples {

// The Detray host store that is used to store the detector
// and the associated memory resource
template <typename memory_source_t>
struct DetrayStore {
  // Constructor from arguments
  DetrayStore(std::shared_ptr<memory_source_t> mSource,
              ActsPlugins::DetrayHostDetector&& det)
      : memoryResource(std::move(mSource)), detector(std::move(det)) {}

  // The memory resource
  std::shared_ptr<memory_source_t> memoryResource = nullptr;
  // The detray detector instance
  ActsPlugins::DetrayHostDetector detector;
};

using DetrayHostStore = DetrayStore<vecmem::host_memory_resource>;

}  // namespace ActsExamples
