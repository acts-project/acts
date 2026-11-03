// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/examples/utils/printable.hpp"

// Boost include(s).
#include <boost/program_options.hpp>

// System include(s).
#include <iosfwd>
#include <memory>
#include <string>
#include <string_view>

namespace traccc::opts {

/// Common base class / interface for all of the program option classes
class interface {
 public:
  /// Constructor on top of a common @c boost::program_options object
  ///
  /// @param description The description of this program option group
  ///
  interface(std::string_view description);
  /// Virtual destructor
  virtual ~interface() = default;

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  virtual void read(const boost::program_options::variables_map& vm);

  /// Helper for turning an interface into a printable set of options.
  virtual std::unique_ptr<configuration_printable> as_printable() const = 0;

  /// Get the description of this program option group
  const boost::program_options::options_description& options() const;

 protected:
  /// (Boost) Description of this program option group
  boost::program_options::options_description m_desc;

  /// (String) Description of this program option group
  std::string m_description;

};  // class interface
}  // namespace traccc::opts
