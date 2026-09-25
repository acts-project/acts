// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Mille/detail/WrappedFileHandle.hpp"

#include <cstdio>
#include <filesystem>
#include <utility>

namespace ActsPlugins {

/// @brief helper to wrap a file handle
WrappedFileHandle::WrappedFileHandle(const std::filesystem::path& outf)
    : m_path(outf) {
  if (!outf.empty()) {
    m_handle = std::fopen(outf.c_str(), "w");
  }
}

WrappedFileHandle::~WrappedFileHandle() {
  if (m_handle != nullptr) {
    std::fclose(m_handle);
  }
}

WrappedFileHandle::WrappedFileHandle(WrappedFileHandle&& other) noexcept
    : m_handle(std::exchange(other.m_handle, nullptr)),
      m_path(std::move(other.m_path)) {}

WrappedFileHandle& WrappedFileHandle::operator=(
    WrappedFileHandle&& other) noexcept {
  if (this != &other) {
    if (m_handle != nullptr) {
      std::fclose(m_handle);
    }
    m_handle = std::exchange(other.m_handle, nullptr);
    m_path = std::move(other.m_path);
  }

  return *this;
}
FILE* WrappedFileHandle::operator()() const {
  return m_handle;
}
bool WrappedFileHandle::isRedirected() const {
  return m_handle != nullptr;
}
const std::filesystem::path& WrappedFileHandle::path() const {
  return m_path;
}

}  // namespace ActsPlugins
