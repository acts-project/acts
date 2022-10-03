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

#include <optional>
#include <string>
#include <vector>

namespace Acts {

struct ExaTrkXTime;

/// @brief Base class for all implementations of the Exa.TrkX pipeline
///
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
  /// @param inputValues Spacepoint data as a flattened NxD array, where D is
  /// the dimensionality of a spacepoint (usually 3, but additional information
  /// like cell information can be provided).
  /// @param spacepointIDs The corresponding spacepoint IDs
  /// @param trackCandidates This vector is filled with the tracks as vectors
  /// of spacepoint IDs
  /// @param logger If provided, logging is enabled
  /// @param recordTiming If enabled, returns a @ref ExaTrkXTime object with
  /// measured timings
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
