// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/ContextualDetector/PrealignedDecorator.hpp"

#include <random>
#include <thread>

ActsExamples::Contextual::PrealignedDecorator::PrealignedDecorator(
    const ActsExamples::Contextual::PrealignedDecorator::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

ActsExamples::ProcessCode
ActsExamples::Contextual::PrealignedDecorator::decorate(
    AlgorithmContext& context) {
  ACTS_VERBOSE("IOV handling in thread " << std::this_thread::get_id() << ".");

  // Create an algorithm local random number generator
  RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);
  std::uniform_int_distribution<unsigned int> iovGenerator(0, m_cfg.nIovs);

  unsigned int iov = iovGenerator(rng);

  // Set the geometry context
  PrealignedDetectorElement::ContextType alignedContext{iov};
  context.geoContext = alignedContext;

  return ProcessCode::SUCCESS;
}
