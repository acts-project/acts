// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepGeometryContext.hpp"

const Acts::Transform3& Acts::DD4hepGeometryContext::contextualTransform(
    const Acts::DD4hepDetectorElement& dElement) const {
  return dElement.transform(DD4hepGeometryContext::inactive());
}
