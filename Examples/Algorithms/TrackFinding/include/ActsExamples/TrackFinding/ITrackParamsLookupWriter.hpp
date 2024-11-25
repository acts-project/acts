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

/// @brief Interface for writing track parameter lookup tables
class ITrackParamsLookupWriter {
 public:
  /// Virtual Destructor
  virtual ~ITrackParamsLookupWriter() = default;

  /// Writer method
  ///
  /// @param lookup track lookup to write
  virtual void writeLookup(const TrackParamsLookup& lookup) const = 0;
};

}  // namespace ActsExamples
