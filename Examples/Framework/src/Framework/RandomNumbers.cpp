// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//
//  RandomNumbers.cpp
//  ACTFW
//
//  Created by Andreas Salzburger on 17/05/16.
//
//

#include "ACTFW/Framework/RandomNumbers.hpp"

FW::RandomNumbers::RandomNumbers(const Config& cfg) : m_cfg(cfg) {}

FW::RandomEngine FW::RandomNumbers::spawnGenerator(
    const AlgorithmContext& context) const {
  return RandomEngine(generateSeed(context));
}

uint64_t FW::RandomNumbers::generateSeed(
    const AlgorithmContext& context) const {
  // use Cantor pairing function to generate a unique generator id from
  // algorithm and event number to get a consistent seed
  // see https://en.wikipedia.org/wiki/Pairing_function#Cantor_pairing_function
  const uint64_t k1 = context.algorithmNumber;
  const uint64_t k2 = context.eventNumber;
  const uint64_t id = (k1 + k2) * (k1 + k2 + 1) / 2 + k2;
  return m_cfg.seed + id;
}
