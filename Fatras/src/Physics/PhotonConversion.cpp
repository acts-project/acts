// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/ElectroMagnetic/PhotonConversion.hpp"

#include "Acts/Definitions/ParticleData.hpp"

const ActsFatras::PhotonConversion::Scalar
    ActsFatras::PhotonConversion::kElectronMass =
        Acts::findMass(Acts::PdgParticle::eElectron).value();
