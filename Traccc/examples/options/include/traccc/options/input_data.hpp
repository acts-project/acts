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
#include <cstddef>
#include <string>

namespace traccc::opts {

/// Options for the input data that a given application would use
class input_data : public interface {
 public:
  /// @name Options
  /// @{

  // Use acts geometry source
  bool use_acts_geom_source = true;
  /// The data format of the input files
  traccc::data_format format = data_format::csv;
  /// Directory of the input files
  std::string directory = "odd/geant4_10muon_10GeV/";
  /// The number of events to process
  std::size_t events = 1;
  /// The number of events to skip
  std::size_t skip = 0;

  /// @}

  /// Constructor
  input_data();

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // class input_data

}  // namespace traccc::opts
