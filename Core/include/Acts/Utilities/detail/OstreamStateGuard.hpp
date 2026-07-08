// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <ios>
#include <ostream>

namespace Acts::detail {

/// @brief Restore an output stream's formatting state on scope exit.
///
/// This keeps diagnostic output helpers from leaking persistent manipulators
/// such as std::fixed, std::setprecision, std::setfill, or std::boolalpha into
/// the stream provided by the caller.
class OstreamStateGuard {
 public:
  explicit OstreamStateGuard(std::ostream& stream)
      : m_stream(&stream),
        m_flags(stream.flags()),
        m_precision(stream.precision()),
        m_width(stream.width()),
        m_fill(stream.fill()) {}

  ~OstreamStateGuard() noexcept {
    m_stream->flags(m_flags);
    m_stream->precision(m_precision);
    m_stream->width(m_width);
    m_stream->fill(m_fill);
  }

  OstreamStateGuard(const OstreamStateGuard&) = delete;
  OstreamStateGuard& operator=(const OstreamStateGuard&) = delete;
  OstreamStateGuard(OstreamStateGuard&&) = delete;
  OstreamStateGuard& operator=(OstreamStateGuard&&) = delete;

 private:
  std::ostream* m_stream{};
  std::ios_base::fmtflags m_flags{};
  std::streamsize m_precision{};
  std::streamsize m_width{};
  char m_fill{};
};

}  // namespace Acts::detail
