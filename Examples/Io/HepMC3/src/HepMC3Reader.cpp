// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/HepMC3/HepMC3Reader.hpp"

bool FW::HepMC3ReaderAscii::readEvent(HepMC3::ReaderAscii& reader,
                                      std::shared_ptr<HepMC3::GenEvent> event) {
  // Read event and store it
  return reader.read_event(*event);
}

bool FW::HepMC3ReaderAscii::status(HepMC3::ReaderAscii& reader) {
  return !reader.failed();
}
