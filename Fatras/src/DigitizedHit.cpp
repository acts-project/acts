// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/DigitizedHit.hpp"

ActsFatras::DigitizedHit::DigitizedHit(
    std::vector<Hit>&& hits, const Acts::Surface& surface,
    std::unique_ptr<const DigitizedHit::IContent> content)

    : m_simulatedHits(std::move(hits)),
      m_surface(surface.getSharedPtr()),
      m_content(std::move(content)) {}

ActsFatras::DigitizedHit::DigitizedHit(const DigitizedHit& other)
    : m_simulatedHits(other.m_simulatedHits),
      m_surface(other.m_surface),
      m_content(nullptr) {}

bool ActsFatras::DigitizedHit::operator==(const DigitizedHit& other) const {
  if (&other != this) {
    return (m_simulatedHits == other.m_simulatedHits and m_surface == other.m_surface);
  }
  return true;
}
