// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

/// Forward declaration of incomplete type cl::sycl::queue
inline namespace cl {
namespace sycl {
class queue;
}
};  // namespace cl

namespace Acts::Sycl {
/// @brief This function creates the SYCL queue object.
///
/// SYCL implementation details are hidden from this class, that's why we're
/// creating the queue in a standalone translation unit.
///
/// @return A pointer to the queue.
std::shared_ptr<cl::sycl::queue> createQueue(const std::string&);
}  // namespace Acts::Sycl
