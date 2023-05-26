// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/StackTrace.hpp"

#include <memory>

#include <boost/stacktrace/detail/frame_unwind.ipp>  // needed to avoid linker error
#include <boost/stacktrace/stacktrace.hpp>

namespace Acts {

struct StackTrace::Impl {
  Impl(boost::stacktrace::stacktrace st) : m_st(std::move(st)) {}
  boost::stacktrace::stacktrace m_st;
};

StackTrace::StackTrace(std::size_t skip, std::size_t maxDepth) {
  m_impl =
      std::make_unique<Impl>(boost::stacktrace::stacktrace{skip + 1, maxDepth});
}

StackTrace::StackTrace(StackTrace &&other) = default;
StackTrace &StackTrace::operator=(StackTrace &&other) = default;

StackTrace::StackTrace(const StackTrace &other) {
  m_impl = std::make_unique<Impl>(*other.m_impl);
}

StackTrace &StackTrace::operator=(const StackTrace &other) {
  m_impl = std::make_unique<Impl>(*other.m_impl);
  return *this;
}

StackTrace::~StackTrace() = default;

std::ostream &operator<<(std::ostream &os, const StackTrace &st) {
  os << st.m_impl->m_st;
  return os;
}

std::pair<std::string, std::size_t> StackTrace::topSourceLocation() const {
  const auto &frame = m_impl->m_st.as_vector().at(0);
  return {frame.source_file(), frame.source_line()};
}

bool operator==(const StackTrace &lhs, const StackTrace &rhs) {
  const auto &fl = *lhs.m_impl->m_st.begin();
  const auto &fr = *rhs.m_impl->m_st.begin();

  return fl.source_file() == fr.source_file() &&
         fl.source_line() == fr.source_line();
}

}  // namespace Acts
