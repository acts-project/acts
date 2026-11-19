// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <detray/core/detector_buffer.hpp>

#include "traccc/geometry/detector_type_list.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/geometry/move_only_any.hpp"

namespace traccc {
using detray::detector_buffer;
using detray::detector_buffer_visitor;

/// Buffer creation follows the detector types enabled for traccc's backends.
inline detector_buffer buffer_from_host_detector(const host_detector& detector,
                                                 vecmem::memory_resource& mr,
                                                 vecmem::copy& copy) {
  return detray::buffer_from_host_detector<detector_type_list>(detector, mr,
                                                               copy);
}
}  // namespace traccc
