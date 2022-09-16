// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/ExaTrkX/ExaTrkXTiming.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <string>
#include <vector>

namespace Acts {

struct ExaTrkXTime;

/// @class ExaTrkXTrackFindingBase
///
/// @brief Base class for all implementations of the Exa.TrkX pipeline
class ExaTrkXTrackFindingBase {
 public:
  /// Constructor
  ///
  /// @param name An identifier available e.g. for debug output
  ExaTrkXTrackFindingBase(const std::string& name) : m_name(name) {}

  /// Destructor
  virtual ~ExaTrkXTrackFindingBase() {}

  ExaTrkXTrackFindingBase() = delete;
  ExaTrkXTrackFindingBase(const ExaTrkXTrackFindingBase&) = delete;
  ExaTrkXTrackFindingBase& operator=(const ExaTrkXTrackFindingBase&) = delete;

  /// Run the inference
  ///
  /// @param inputValues tPacked spacepoints in the form
  /// [ r1, phi1, z1, r2, phi2, z2, ... ]
  /// @param spacepointIDs The corresponding spacepoint spacepoint spacepointIDs
  /// @param trackCandidates This vector is filled with the tracks as vectors of spacepoint spacepoint IDs
  /// @param logger If provided, logging is enabled
  /// @note The input values are not const, because the ONNX API
  /// takes only non-const pointers.
  virtual std::optional<ExaTrkXTime> getTracks(
      std::vector<float>& inputValues, std::vector<int>& spacepointIDs,
      std::vector<std::vector<int> >& trackCandidates,
      Acts::LoggerWrapper logger = Acts::getDummyLogger(),
      bool recordTiming = false) const = 0;

  /// Returns the name of the algorithm
  const std::string& name() const { return m_name; }

 private:
  std::string m_name;
};

}  // namespace Acts
