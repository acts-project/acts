// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsProfiler.hpp"

#include <algorithm>
#include <memory>
#include <utility>

namespace Acts::Experimental {

void GbtsProfiler::setEventNumber(const std::uint64_t eventNumber) {
  const std::scoped_lock lock(m_mutex);
  m_eventNumbers[std::this_thread::get_id()] = eventNumber;
}

void GbtsProfiler::add(Event&& event) {
  const std::scoped_lock lock(m_mutex);
  event.eventNumber = m_eventNumbers[std::this_thread::get_id()];
  m_events.push_back(std::move(event));
}

std::unique_ptr<const GbtsProfiler::Event> GbtsProfiler::take(
    const std::uint64_t eventNumber) {
  const std::scoped_lock lock(m_mutex);
  const auto it = std::ranges::find(m_events, eventNumber, &Event::eventNumber);
  if (it == m_events.end()) {
    return nullptr;
  }
  auto event = std::make_unique<const Event>(std::move(*it));
  m_events.erase(it);
  return event;
}

std::vector<std::unique_ptr<const GbtsProfiler::Event>>
GbtsProfiler::takeAll() {
  const std::scoped_lock lock(m_mutex);
  std::vector<std::unique_ptr<const Event>> events;
  events.reserve(m_events.size());
  for (Event& event : m_events) {
    events.push_back(std::make_unique<const Event>(std::move(event)));
  }
  m_events.clear();
  return events;
}

}  // namespace Acts::Experimental
