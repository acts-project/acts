// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <utility>

namespace ActsExamples {

/// Helper that keeps ROOT branch payload alive via unique_ptr while also
/// providing the raw pointer/pointer-to-pointer interface ROOT expects.
template <typename T>
class RootBranchPtr {
 public:
  using value_type = T;

  RootBranchPtr() = default;
  explicit RootBranchPtr(std::nullptr_t) {}

  RootBranchPtr(const RootBranchPtr&) = delete;
  RootBranchPtr& operator=(const RootBranchPtr&) = delete;

  RootBranchPtr(RootBranchPtr&& other) noexcept { *this = std::move(other); }

  RootBranchPtr& operator=(RootBranchPtr&& other) noexcept {
    if (this != &other) {
      m_storage = std::move(other.m_storage);
      m_pointer = m_storage ? m_storage.get() : nullptr;
      other.m_pointer = nullptr;
    }
    return *this;
  }

  /// Allocate owned storage if not already available.
  void allocate() {
    if (!m_storage) {
      m_storage = std::make_unique<T>();
      m_pointer = m_storage.get();
    }
  }

  /// Release owned storage and clear the raw pointer.
  void reset() {
    m_storage.reset();
    m_pointer = nullptr;
  }

  bool hasValue() const { return m_pointer != nullptr; }

  T& operator*() { return *m_pointer; }
  const T& operator*() const { return *m_pointer; }

  T* operator->() { return m_pointer; }
  const T* operator->() const { return m_pointer; }

  /// Returns a reference to the managed pointer so callers can take its
  /// address when interfacing with ROOT.
  T*& get() { return m_pointer; }
  T* const& get() const { return m_pointer; }

 private:
  std::unique_ptr<T> m_storage = nullptr;
  T* m_pointer = nullptr;
};

}  // namespace ActsExamples
