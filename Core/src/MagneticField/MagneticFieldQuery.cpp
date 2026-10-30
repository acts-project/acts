// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldQuery.hpp"

#include <stdexcept>
#include <utility>

Acts::MagneticFieldQuery::MagneticFieldQuery(
    std::shared_ptr<const MagneticFieldProvider> provider,
    const MagneticFieldContext& mctx)
    : m_provider(std::move(provider)) {
  if (!m_provider) {
    throw std::invalid_argument("MagneticFieldQuery: provider is null");
  }
  m_cache = m_provider->makeCache(mctx);
}

Acts::Result<Acts::Vector3> Acts::MagneticFieldQuery::getField(
    const Vector3& position) {
  return m_provider->getField(position, m_cache);
}
