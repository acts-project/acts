// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include <cstdio>
#include <filesystem>

namespace ActsPlugins::ActsToMille {

/// @brief helper to wrap a file handle
class WrappedFileHandle {
 public:
  explicit WrappedFileHandle(const std::filesystem::path& outf = "");
  ~WrappedFileHandle();
  WrappedFileHandle(const WrappedFileHandle&) = delete;
  WrappedFileHandle& operator=(const WrappedFileHandle&) = delete;

  WrappedFileHandle(WrappedFileHandle&& other) noexcept;

  WrappedFileHandle& operator=(WrappedFileHandle&& other) noexcept;
  FILE* operator()() const;
  bool isRedirected() const;
  const std::filesystem::path& path() const;

 private:
  FILE* m_handle = nullptr;
  std::filesystem::path m_path{};
};

}  // namespace ActsPlugins::ActsToMille
