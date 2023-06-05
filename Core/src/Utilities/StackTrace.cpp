// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/StackTrace.hpp"

#include <memory>
#include <memory_resource>

#include <boost/stacktrace/detail/frame_unwind.ipp>  // needed to avoid linker error
#include <boost/stacktrace/frame.hpp>
#include <boost/stacktrace/stacktrace.hpp>

namespace {
using stacktrace_t = boost::stacktrace::basic_stacktrace<
    std::pmr::polymorphic_allocator<boost::stacktrace::frame>>;
}  // namespace

namespace Acts {

struct StackTrace::Impl {
  Impl(stacktrace_t st) : m_st(std::move(st)) {}
  stacktrace_t m_st;
};

StackTrace::StackTrace(std::size_t skip, std::size_t maxDepth,
                       std::pmr::memory_resource &mem)
    : m_allocator{&mem} {
  m_impl = m_allocator.allocate(1);
  m_allocator.construct(m_impl, stacktrace_t{skip + 1, maxDepth, &mem});
}

StackTrace::StackTrace(const StackTrace &other)
    : m_allocator{other.m_allocator.resource()} {
  m_impl = m_allocator.allocate(1);
  m_allocator.construct(m_impl, *other.m_impl);
}

StackTrace &StackTrace::operator=(const StackTrace &other) {
  m_impl = m_allocator.allocate(1);
  m_allocator.construct(m_impl, *other.m_impl);
  return *this;
}

StackTrace::~StackTrace() {
  assert(m_impl != nullptr);
  m_allocator.destroy(m_impl);
}

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

  // @TODO: Check all frames
  return boost::stacktrace::hash_value(fl) == boost::stacktrace::hash_value(fr);
}

}  // namespace Acts
