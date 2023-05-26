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
#include <string>

namespace Acts {
struct StackTrace {
  StackTrace(std::size_t skip, std::size_t maxDepth);
  StackTrace(StackTrace &&other);
  StackTrace(const StackTrace &other);
  StackTrace &operator=(StackTrace &&other);
  StackTrace &operator=(const StackTrace &other);

  std::pair<std::string, std::size_t> topSourceLocation() const;

  ~StackTrace();

  friend std::ostream &operator<<(std::ostream &os, const StackTrace &st);

  friend bool operator==(const StackTrace &lhs, const StackTrace &rhs);

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace Acts
