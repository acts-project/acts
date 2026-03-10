// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryModule.h"

#include <filesystem>
#include <memory>

namespace Acts {

/// Returns the host ABI tag expected by the runtime module loader.
const char* geometryModuleHostAbiTag() noexcept;

/// Opaque runtime module build handle.
class GeometryModuleHandle {
 public:
  GeometryModuleHandle() = default;
  ~GeometryModuleHandle();

  GeometryModuleHandle(const GeometryModuleHandle&) = delete;
  GeometryModuleHandle& operator=(const GeometryModuleHandle&) = delete;
  GeometryModuleHandle(GeometryModuleHandle&& other) noexcept;
  GeometryModuleHandle& operator=(GeometryModuleHandle&& other) noexcept;

  void* get() const noexcept { return m_handle; }
  explicit operator bool() const noexcept { return m_handle != nullptr; }

 private:
  using DestroyFn = void (*)(void*);

  GeometryModuleHandle(void* handle, DestroyFn destroy,
                       std::shared_ptr<void> library) noexcept;
  void reset() noexcept;

  void* m_handle{nullptr};
  DestroyFn m_destroy{nullptr};
  std::shared_ptr<void> m_library;

  friend GeometryModuleHandle loadGeometryModule(
      const std::filesystem::path& modulePath);
};

/// Load a module shared library, validate ABI compatibility, build and return
/// an opaque handle.
GeometryModuleHandle loadGeometryModule(
    const std::filesystem::path& modulePath);

}  // namespace Acts
