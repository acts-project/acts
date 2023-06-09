// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/PhysicsListFactory.hpp"

#include <stdexcept>

#include <FTFP_BERT.hh>
#include <FTFP_BERT_ATL.hh>

namespace ActsExamples {

G4VModularPhysicsList* PhysicsListFactory::factorize(
    const std::string& list) const {
  if (list == "FTFP_BERT") {
    return new FTFP_BERT();
  } else if (list == "FTFP_BERT_ATL") {
    return new FTFP_BERT_ATL();
  }

  throw std::invalid_argument("unhandled list");
}

}  // namespace ActsExamples
