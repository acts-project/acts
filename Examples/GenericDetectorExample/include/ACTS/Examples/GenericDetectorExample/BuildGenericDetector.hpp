// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_BUILDGENERICDETECTOR_H
#define ACTS_BUILDGENERICDETECTOR_H 1

// STL include(s)
#include <memory>

// ACTS include(s)
namespace Acts
{
  class TrackingGeometry;

  std::unique_ptr<const Acts::TrackingGeometry> trackingGeometry();
}

#endif // ACTS_BUILDGENERICDETECTOR_H
