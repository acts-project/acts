// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/io/data_format.hpp"
#include "traccc/options/details/interface.hpp"

// System include(s).
#include <string>
#include <string_view>

namespace traccc::opts {

/// Options for the output data that a given application would produce
class output_data : public interface {
 public:
  /// @name Options
  /// @{

  /// The data format of the input files
  traccc::data_format format = data_format::csv;
  /// Directory of the input files
  std::string directory = "testing/";

  /// @}

  /// Constructor, with default arguments
  ///
  /// @param format The data format for the output
  /// @param directory The directory for the output
  ///
  output_data(traccc::data_format format = data_format::csv,
              std::string_view directory = "testing/");

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // class output_data

}  // namespace traccc::opts
