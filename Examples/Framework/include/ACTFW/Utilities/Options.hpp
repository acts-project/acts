// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iosfwd>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace FW {
namespace Options {

/// Half open [lower,upper) interval type for user options.
///
/// A missing limit represents an unbounded upper or lower limit. With just
/// one defined limit the interval is just a lower/upper bound; with both
/// limits undefined, the interval is unbounded everywhere and thus contains
/// all possible values.
///
/// This is intended as a utility type for the user options and not as a
/// variable type for the configuration structs. Simple primitive types should
/// be preferred there.
struct Interval {
  std::optional<double> lower;
  std::optional<double> upper;
};

/// Extract an interval from an input of the form 'lower:upper'.
///
/// An input of the form `lower:` or `:upper` sets just one of the limits. Any
/// other input leads to an unbounded interval. If the input is `:SECOND` the
///
/// @note The more common range notation uses `lower-upper` but the `-`
///       separator complicates the parsing of negative values.
std::istream& operator>>(std::istream& is, Interval& interval);

/// Print an interval as `lower:upper`.
std::ostream& operator<<(std::ostream& os, const Interval& interval);

}  // namespace Options
}  // namespace FW

using read_series = std::vector<int>;
using read_range = std::vector<double>;
using read_strings = std::vector<std::string>;

// Overloads must exist in the `std` namespace so ADL-lookup can find them.
namespace std {

std::ostream& operator<<(std::ostream& os, const read_series& vec);

std::ostream& operator<<(std::ostream& os, const read_range& vec);

std::ostream& operator<<(std::ostream& os, const read_strings& vec);

}  // namespace std
