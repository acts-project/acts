// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/RandomNumbers.hpp"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <boost/functional/hash.hpp>

namespace ActsExamples {

RandomNumberChannel::RandomNumberChannel(const RandomNumbers& randomNumbers,
                                         RandomSeed seed)
    : m_randomNumbers(&randomNumbers), m_seed(seed) {}

RandomSeed RandomNumberChannel::seed() const {
  return m_seed;
}

RandomEngine RandomNumberChannel::createEngine() const {
  return RandomEngine(m_seed);
}

RandomNumberChannel RandomNumberChannel::createSubChannel(
    RandomSeed seed) const {
  return RandomNumberChannel(*m_randomNumbers, m_seed + seed);
}

RandomNumbers::RandomNumbers(const Config& cfg) : m_cfg(cfg) {}

RandomNumberChannel RandomNumbers::createChannel() const {
  return RandomNumberChannel(*this, m_cfg.seed);
}

RandomNumberChannel RandomNumbers::createAlgorithmEventChannel(
    const AlgorithmContext& context) const {
  return RandomNumberChannel(*this, generateSeed(context));
}

RandomEngine RandomNumbers::spawnGenerator(
    const AlgorithmContext& context) const {
  return RandomEngine(generateSeed(context));
}

RandomSeed RandomNumbers::generateSeed(const AlgorithmContext& context) const {
  return m_cfg.seed + context.eventNumber;
}

}  // namespace ActsExamples
