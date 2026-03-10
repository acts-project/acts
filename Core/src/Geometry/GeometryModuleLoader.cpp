// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleLoader.hpp"

#if defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h>
#endif

#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

#ifndef ACTS_GEOMETRY_MODULE_ABI_TAG
#error "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided by CMake when building ActsCore."
#endif

namespace {

#if defined(__unix__) || defined(__APPLE__)
using GeometryModuleEntryPointV1 = const ActsGeometryModuleV1* (*)(void);

std::shared_ptr<void> openSharedLibrary(const std::filesystem::path& path) {
  void* rawHandle = ::dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (rawHandle == nullptr) {
    const char* error = ::dlerror();
    std::ostringstream msg;
    msg << "Failed to load geometry module '" << path.string() << "'";
    if (error != nullptr) {
      msg << ": " << error;
    }
    throw std::runtime_error(msg.str());
  }
  return std::shared_ptr<void>(rawHandle, [](void* handle) {
    if (handle != nullptr) {
      ::dlclose(handle);
    }
  });
}

GeometryModuleEntryPointV1 resolveEntrypointV1(
    const std::filesystem::path& path, const std::shared_ptr<void>& library) {
  ::dlerror();
  void* symbol = ::dlsym(library.get(), "acts_geometry_module_v1");
  const char* error = ::dlerror();
  if (error != nullptr) {
    std::ostringstream msg;
    msg << "Failed to resolve acts_geometry_module_v1 in '" << path.string()
        << "': " << error;
    throw std::runtime_error(msg.str());
  }
  if (symbol == nullptr) {
    std::ostringstream msg;
    msg << "Entry point acts_geometry_module_v1 resolved to nullptr in '"
        << path.string() << "'";
    throw std::runtime_error(msg.str());
  }
  return reinterpret_cast<GeometryModuleEntryPointV1>(symbol);
}

#endif

}  // namespace

namespace Acts {

const char* geometryModuleHostAbiTag() noexcept {
  return ACTS_GEOMETRY_MODULE_ABI_TAG;
}

ActsGeometryModuleRequestV1 makeGeometryModuleRequest(void* context,
                                                      void* userData) noexcept {
  return ActsGeometryModuleRequestV1{
      .abi_tag = geometryModuleHostAbiTag(),
      .context = context,
      .user_data = userData,
  };
}

GeometryModuleHandle::GeometryModuleHandle(
    void* handle, DestroyFn destroy, std::shared_ptr<void> library) noexcept
    : m_handle(handle), m_destroy(destroy), m_library(std::move(library)) {}

GeometryModuleHandle::~GeometryModuleHandle() { reset(); }

GeometryModuleHandle::GeometryModuleHandle(GeometryModuleHandle&& other) noexcept
    : m_handle(std::exchange(other.m_handle, nullptr)),
      m_destroy(std::exchange(other.m_destroy, nullptr)),
      m_library(std::move(other.m_library)) {}

GeometryModuleHandle& GeometryModuleHandle::operator=(
    GeometryModuleHandle&& other) noexcept {
  if (this == &other) {
    return *this;
  }
  reset();
  m_handle = std::exchange(other.m_handle, nullptr);
  m_destroy = std::exchange(other.m_destroy, nullptr);
  m_library = std::move(other.m_library);
  return *this;
}

void GeometryModuleHandle::reset() noexcept {
  if (m_handle != nullptr && m_destroy != nullptr) {
    m_destroy(m_handle);
  }
  m_handle = nullptr;
  m_destroy = nullptr;
  m_library.reset();
}

GeometryModuleHandle loadGeometryModule(
    const std::filesystem::path& modulePath,
    const ActsGeometryModuleRequestV1& request) {
#if !(defined(__unix__) || defined(__APPLE__))
  (void)modulePath;
  (void)request;
  throw std::runtime_error(
      "Runtime geometry modules are only supported on Unix-like systems");
#else
  if (request.abi_tag == nullptr) {
    throw std::invalid_argument("Geometry module request ABI tag is null");
  }
  if (std::strcmp(request.abi_tag, geometryModuleHostAbiTag()) != 0) {
    std::ostringstream msg;
    msg << "Geometry module host ABI tag mismatch: request='" << request.abi_tag
        << "' host='" << geometryModuleHostAbiTag() << "'";
    throw std::runtime_error(msg.str());
  }

  auto library = openSharedLibrary(modulePath);
  auto entryPoint = resolveEntrypointV1(modulePath, library);
  const ActsGeometryModuleV1* descriptor = entryPoint();
  if (descriptor == nullptr) {
    throw std::runtime_error("Geometry module descriptor is null");
  }
  if (descriptor->module_abi_tag == nullptr || descriptor->build == nullptr ||
      descriptor->destroy == nullptr || descriptor->last_error == nullptr) {
    throw std::runtime_error("Geometry module descriptor is incomplete");
  }
  if (std::strcmp(descriptor->module_abi_tag, request.abi_tag) != 0) {
    std::ostringstream msg;
    msg << "Geometry module ABI mismatch: module='"
        << descriptor->module_abi_tag << "' host='" << request.abi_tag
        << "' path='" << modulePath.string() << "'";
    throw std::runtime_error(msg.str());
  }

  void* rawHandle = descriptor->build(&request);
  if (rawHandle == nullptr) {
    const char* moduleError = descriptor->last_error();
    std::ostringstream msg;
    msg << "Geometry module build returned null handle";
    if (moduleError != nullptr && moduleError[0] != '\0') {
      msg << ": " << moduleError;
    }
    throw std::runtime_error(msg.str());
  }
  return GeometryModuleHandle(rawHandle, descriptor->destroy, library);
#endif
}

}  // namespace Acts
