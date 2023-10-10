// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/FpeMonitoring/FpeMonitor.hpp"
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <memory>

namespace ActsExamples {

class WhiteBoard;

/// Aggregated information to run one algorithm over one event.
struct AlgorithmContext {
  /// @brief constructor with arguments
  ///
  /// @param alg is the algorithm/service/writer number
  /// @param event is the event number
  /// @param store is the event-wise event store
  ///
  /// @note the event dependent contexts are to be added by the
  /// Sequencer::m_decorators list
  AlgorithmContext(size_t alg, size_t event, WhiteBoard& store)
      : algorithmNumber(alg), eventNumber(event), eventStore(store) {}

  /// @brief ++operator overload to increase the algorithm number
  AlgorithmContext& operator++() {
    algorithmNumber += 1;
    return (*this);
  }

  size_t algorithmNumber;            ///< Unique algorithm identifier
  size_t eventNumber;                ///< Unique event identifier
  WhiteBoard& eventStore;            ///< Per-event data store
  Acts::GeometryContext geoContext;  ///< Per-event geometry context
  Acts::MagneticFieldContext
      magFieldContext;                    ///< Per-event magnetic Field context
  Acts::CalibrationContext calibContext;  ///< Per-event calibration context

  Acts::FpeMonitor* fpeMonitor = nullptr;
};

}  // namespace ActsExamples
