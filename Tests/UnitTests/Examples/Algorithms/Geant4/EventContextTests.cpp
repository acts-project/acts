// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Geant4/EventStore.hpp"

#include <stdexcept>
#include <utility>

namespace ActsExamples::Geant4 {

BOOST_AUTO_TEST_CASE(BindingEndsOnReturnAndException) {
  const auto geometry = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const SimParticleContainer particles;
  EventContext context;
  const EventContext& reader = context;
  BOOST_CHECK_THROW(reader.current(), std::logic_error);

  auto transport = [&](bool fail) -> EventResult {
    EventStore event(geometry, particles);
    auto binding = context.bind(event);
    BOOST_CHECK_EQUAL(&reader.current().inputParticles, &particles);
    BOOST_CHECK_THROW((void)context.bind(event), std::logic_error);
    context.current().result.numberGeantSteps = 1;
    if (fail) {
      throw std::runtime_error("transport failed");
    }
    return std::move(event.result);
  };

  BOOST_CHECK_THROW(transport(true), std::runtime_error);
  BOOST_CHECK_THROW(context.current(), std::logic_error);
  auto result = transport(false);
  BOOST_CHECK_THROW(reader.current(), std::logic_error);
  BOOST_CHECK_EQUAL(result.numberGeantSteps, 1u);
}

}  // namespace ActsExamples::Geant4
