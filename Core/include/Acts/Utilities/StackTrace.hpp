// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <memory>
#include <memory_resource>
#include <string>

namespace Acts {
struct StackTrace {
  StackTrace(std::size_t skip, std::size_t maxDepth,
             std::pmr::memory_resource &mem = *std::pmr::new_delete_resource());
  StackTrace(const StackTrace &other);
  StackTrace &operator=(const StackTrace &other);

  std::pair<std::string, std::size_t> topSourceLocation() const;

  ~StackTrace();

  friend std::ostream &operator<<(std::ostream &os, const StackTrace &st);

  friend bool operator==(const StackTrace &lhs, const StackTrace &rhs);

 private:
  struct Impl;
  std::pmr::polymorphic_allocator<Impl> m_allocator;

  Impl *m_impl{nullptr};
};

}  // namespace Acts
