// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#ifdef ACTS_EXATRKX_NVTX
#include <optional>

#include <nvtx3/nvToolsExt.h>

namespace Acts::detail {

struct NvtxRange {
  nvtxRangeId_t id;
  NvtxRange(const char *name) { id = nvtxRangeStartA(name); }
  NvtxRange(const NvtxRange &) = delete;
  NvtxRange(NvtxRange &&) = delete;

  ~NvtxRange() { nvtxRangeEnd(id); }
};

}  // namespace Acts::detail

#define ACTS_NVTX_START(name)                                \
  std::optional<Acts::detail::NvtxRange> _nvtx##name(#name); \
  do {                                                       \
  } while (0)
#define ACTS_NVTX_STOP(name) \
  _nvtx##name.reset();       \
  do {                       \
  } while (0)
#else
#define ACTS_NVTX_START(name) \
  do { /*nothing*/            \
  } while (0)
#define ACTS_NVTX_STOP(name) \
  do { /*nothing*/           \
  } while (0)
#endif
