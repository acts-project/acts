// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/ProtoVolumeMaterial.hpp"

#include <ostream>

Acts::ProtoVolumeMaterial::ProtoVolumeMaterial(const MultiAxisSpec3D& binning)
    : m_binning(binning) {}

std::ostream& Acts::ProtoVolumeMaterial::toStream(std::ostream& stream) const {
  return stream << "Acts::ProtoVolumeMaterial : " << m_binning;
}
