// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/StackTrace.hpp"

#include <iterator>
#include <memory>
#include <memory_resource>
#include <sstream>

#include <boost/stacktrace/detail/frame_unwind.ipp>  // needed to avoid linker error
#include <boost/stacktrace/frame.hpp>
#include <boost/stacktrace/stacktrace.hpp>

namespace {
using stacktrace_t = boost::stacktrace::basic_stacktrace<
    std::pmr::polymorphic_allocator<boost::stacktrace::frame>>;
}  // namespace

namespace Acts {

struct StackTrace::Impl {
  Impl(std::pmr::memory_resource &mem) : frames{&mem} {}
  std::pmr::vector<boost::stacktrace::frame> frames;
};

StackTrace::StackTrace(std::size_t skip, std::size_t maxDepth,
                       std::pmr::memory_resource &mem)
    : m_allocator{&mem} {
  m_impl = m_allocator.allocate(1);
  m_allocator.construct(m_impl, *m_allocator.resource());

  stacktrace_t st{skip + 1, maxDepth, &mem};
  std::copy(st.begin(), st.end(), std::back_inserter(m_impl->frames));
}

StackTrace::StackTrace(const StackTrace &other)
    : m_allocator{std::pmr::new_delete_resource()} {
  m_impl = m_allocator.allocate(1);
  m_allocator.construct(m_impl, *m_allocator.resource());
  std::copy(other.m_impl->frames.begin(), other.m_impl->frames.end(),
            std::back_inserter(m_impl->frames));
}

StackTrace::StackTrace(StackTrace &&other) : m_allocator{other.m_allocator} {
  m_impl = other.m_impl;
  other.m_impl = nullptr;
}

StackTrace &StackTrace::operator=(const StackTrace &other) {
  assert(m_impl != nullptr);
  std::copy(other.m_impl->frames.begin(), other.m_impl->frames.end(),
            std::back_inserter(m_impl->frames));
  return *this;
}

std::string StackTrace::toString(std::size_t depth) const {
  return boost::stacktrace::detail::to_string(
      m_impl->frames.data(), std::min(depth, m_impl->frames.size()));
}

StackTrace::~StackTrace() {
  if (m_impl != nullptr) {
    m_impl->~Impl();
    m_allocator.deallocate(m_impl, 1);
  }
}

std::ostream &operator<<(std::ostream &os, const StackTrace &st) {
  os << st.toString();
  return os;
}

std::pair<std::string, std::size_t> StackTrace::topSourceLocation() const {
  const auto &frame = m_impl->frames.at(0);
  return {frame.source_file(), frame.source_line()};
}

bool operator==(const StackTrace &lhs, const StackTrace &rhs) {
  const auto &fl = lhs.m_impl->frames;
  const auto &fr = rhs.m_impl->frames;

  // only compare first frame
  return (boost::stacktrace::hash_value(fl.front()) ==
          boost::stacktrace::hash_value(fr.front()));

}

}  // namespace Acts
