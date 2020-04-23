// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/HepMC3/HepMC3Writer.hpp"

bool FW::HepMC3WriterAscii::writeEvent(
    HepMC3::WriterAscii& writer, std::shared_ptr<HepMC3::GenEvent> event) {
  // Write event from storage
  writer.write_event(*event);
  return true;
}

bool FW::HepMC3WriterAscii::status(HepMC3::WriterAscii& writer) {
  return writer.failed();
}
