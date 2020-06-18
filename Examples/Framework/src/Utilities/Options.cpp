// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Utilities/Options.hpp"

#include <cstdlib>
#include <istream>
#include <ostream>

namespace {
static constexpr char s_rangeSeparator = ':';
}

std::istream& FW::Options::operator>>(std::istream& is,
                                      FW::Options::Interval& interval) {
  std::string buf;
  is >> buf;

  // default to an unbounded interval
  interval.lower.reset();
  interval.upper.reset();

  // find the limit separator
  auto pos = buf.find_first_of(s_rangeSeparator);
  // no separator -> invalid input -> unbounded interval
  if (pos == std::string::npos) {
    return is;
  }

  // if it exists, parse limit before separator
  if (0 < pos) {
    auto lowerStr = buf.substr(0, pos);
    interval.lower = std::atof(lowerStr.c_str());
  }
  // if it exists, parse limit after separator
  if ((pos + 1) < buf.size()) {
    auto upperStr = buf.substr(pos + 1);
    interval.upper = std::atof(upperStr.c_str());
  }

  return is;
}

std::ostream& FW::Options::operator<<(std::ostream& os,
                                      const FW::Options::Interval& interval) {
  if (not interval.lower.has_value() and not interval.upper.has_value()) {
    os << "unbounded";
  } else {
    if (interval.lower.has_value()) {
      os << interval.lower.value();
    }
    os << s_rangeSeparator;
    if (interval.upper.has_value()) {
      os << interval.upper.value();
    }
  }
  return os;
}

std::istream& FW::Options::operator>>(
    std::istream& is, std::vector<FW::Options::Interval>& intervals) {
  for (auto& interval : intervals) {
    is >> interval;
  }
  return is;
}

std::ostream& FW::Options::operator<<(
    std::ostream& os, const std::vector<FW::Options::Interval>& intervals) {
  for (auto& interval : intervals) {
    os << interval;
  }
  return os;
}

namespace {
/// Helper function to print multiple elements in a container.
template <typename Iterator>
inline std::ostream& printContainer(std::ostream& os, Iterator begin,
                                    Iterator end, const char* separator) {
  for (auto it = begin; it != end; ++it) {
    if (it != begin) {
      os << separator;
    }
    os << *it;
  }
  return os;
}
}  // namespace

std::ostream& std::operator<<(std::ostream& os, const read_series& vec) {
  return printContainer(os, vec.begin(), vec.end(), " ");
}

std::ostream& std::operator<<(std::ostream& os, const read_range& vec) {
  return printContainer(os, vec.begin(), vec.end(), " ");
}

std::ostream& std::operator<<(std::ostream& os, const read_strings& vec) {
  return printContainer(os, vec.begin(), vec.end(), " ");
}
