// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/FpeMonitoring/MemoryResourceCompat.hpp"

#include <atomic>

#if defined(ACTS_HAVE_PMR_MEMORY_RESOURCE)
namespace std::pmr {
#else
namespace std::experimental {
inline namespace fundamentals_v1 {
namespace pmr {
#endif

namespace {

std::atomic<memory_resource*>& default_resource() noexcept {
  static std::atomic<memory_resource*> res{new_delete_resource()};
  return res;
}

}  // namespace

memory_resource* new_delete_resource() noexcept {
  // vecmem::host_memory_resource is not singleton, and does not check
  // equality via pointer equality, but is close enough for our purposes.
  static vecmem::host_memory_resource res{};
  return &res;
}

memory_resource* get_default_resource() noexcept {
  return default_resource();
}

memory_resource* set_default_resource(memory_resource* res) noexcept {
  memory_resource* new_res = res == nullptr ? new_delete_resource() : res;

  return default_resource().exchange(new_res);
}
#if defined(ACTS_HAVE_PMR_MEMORY_RESOURCE)
}
#else
}  // namespace pmr
}  // namespace fundamentals_v1
}  // namespace std::experimental
#endif
