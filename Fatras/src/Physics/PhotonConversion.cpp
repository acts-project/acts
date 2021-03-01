// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Physics/PhotonConversion/PhotonConversion.hpp"

#include "ActsFatras/Utilities/ParticleData.hpp"

const ActsFatras::PhotonConversion::Scalar
    ActsFatras::PhotonConversion::kElectronMass =
        ActsFatras::findMass(Acts::PdgParticle::eElectron);
