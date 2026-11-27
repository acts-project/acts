// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

#include <traccc/bfield/magnetic_field.hpp>
#include <traccc/utils/memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

namespace ActsExamples {

/// Interface for a traccc device backend (CUDA, HIP, SYCL, Alpaka, ...)
///
/// A backend owns the device memory resources, the copy object and the
/// stream/queue. The code using the backend only ever sees this interface,
/// so it does not depend on any backend headers.
///
/// @note All algorithms share a single stream/queue, so the sequencer has to
///       run single threaded.
class TracccBackend {
 public:
  virtual ~TracccBackend() = default;

  /// Human readable name of the backend
  virtual std::string name() const = 0;

  /// Device (main) and host-accessible memory resources
  virtual const traccc::memory_resource& memoryResource() const = 0;

  /// Copy object for host <-> device transfers on the backend stream/queue
  virtual vecmem::copy& copy() const = 0;

  /// Block until all work submitted to the backend stream/queue is done
  virtual void synchronize() const = 0;

  /// Make a copy of a host magnetic field usable on the device
  virtual traccc::magnetic_field makeMagneticField(
      const traccc::magnetic_field& hostField) const = 0;
};

}  // namespace ActsExamples
