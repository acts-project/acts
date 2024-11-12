// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFinding/TrackParamsLookupTable.hpp"

namespace ActsExamples {

/// @brief Interface for reading track parameter lookup tables
class ITrackParamsLookupReader {
 public:
  /// Virtual Destructor
  virtual ~ITrackParamsLookupReader() = default;

  /// Reader method
  ///
  /// @param path the path to the file to read
  virtual TrackParamsLookup readLookup(const std::string& path) const = 0;
};

}  // namespace ActsExamples
