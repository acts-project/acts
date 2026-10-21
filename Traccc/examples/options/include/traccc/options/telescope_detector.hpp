// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/options/details/interface.hpp"
#include "traccc/options/details/value_array.hpp"

namespace traccc::opts {

/// Command line options used in the telescope detector tests
class telescope_detector : public interface {
 public:
  /// @name Options
  /// @{

  /// Build detector without materials
  bool empty_material = false;
  /// Number of planes
  unsigned int n_planes = 9;
  /// Slab thickness in [mm]
  float thickness = 0.5f;
  /// Space between planes in [mm]
  float spacing = 20.f;
  /// Measurement smearing in [um]
  float smearing = 50.f;
  /// Half length of plane [mm]
  float half_length = 1000000.f;
  /// Vector for plane alignment
  opts::value_array<float, 3> align_vector{0.f, 0.f, 1.f};

  /// @}

  /// Constructor
  telescope_detector();

  /// Read/process the command line options
  ///
  /// @param vm The command line options to interpret/read
  ///
  void read(const boost::program_options::variables_map& vm) override;

  std::unique_ptr<configuration_printable> as_printable() const override;
};  // ckass telescope_detector

}  // namespace traccc::opts
