// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <HepMC3/GenEvent.h>
#include <HepMC3/ReaderAscii.h>

namespace FW {

/// HepMC3 event reader.
struct HepMC3ReaderAscii {
 public:
  /// @brief Reads an event from file
  /// @param reader reader of run files
  /// @param event storage of the read event
  /// @return boolean indicator if the reading was successful
  bool readEvent(HepMC3::ReaderAscii& reader,
                 std::shared_ptr<HepMC3::GenEvent> event);

  /// @brief Reports the status of the reader
  /// @param reader reader of run files
  /// @return boolean status indicator
  bool status(HepMC3::ReaderAscii& reader);
};
}  // namespace FW
