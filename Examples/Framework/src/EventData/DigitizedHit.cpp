// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/DigitizedHit.hpp"

ActsExamples::DigitizedHit::DigitizedHit(const Acts::Surface& surface)
    : m_surface(surface.getSharedPtr()) {}

bool ActsExamples::DigitizedHit::operator==(const DigitizedHit& other) const {
  return (this == &other);
}

const Acts::Surface& ActsExamples::DigitizedHit::referenceSurface() const {
  return (*m_surface.get());
}
