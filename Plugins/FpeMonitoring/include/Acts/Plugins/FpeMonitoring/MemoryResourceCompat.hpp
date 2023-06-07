// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// See https://github.com/acts-project/vecmem/blob/main/core/include/vecmem/memory/memory_resource.hpp

// The purpose of this file is to provide uniform access (on a source-code
// level) to the memory_resource type from the standard library. These are
// either in the std::pmr namespace or in the std::experimental::pmr namespace
// depending on the GCC version used, so we try to unify them by aliassing
// depending on the compiler feature flags.

#if defined(VECMEM_HAVE_PMR_MEMORY_RESOURCE)
#include <memory_resource>

namespace vecmem {
  using memory_resource = std::pmr::memory_resource;
  
}
#elif defined(VECMEM_HAVE_EXPERIMENTAL_PMR_MEMORY_RESOURCE)
#include <experimental/memory_resource>

namespace vecmem {
  using memory_resource = std::experimental::pmr::memory_resource;
  
}
#else
#error "C++17 LFTS V1 (P0220R1) component memory_resource not found!?!"
#endif
