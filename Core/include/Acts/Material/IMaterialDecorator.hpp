// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IMaterialDecorator.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

namespace Acts {

class Surface;
class TrackingVolume;

/// @class IMaterialDecorator
///
/// Virtual base class of surface based material description
//
/// Material associated with a Volume (homogenous, binned, interpolated)
class IMaterialDecorator
{
public:
  /// Virtual Destructor
  virtual ~IMaterialDecorator() = default;

  /// Decorate a surface
  virtual void
  decorate(Surface& surface) const = 0;

  /// Decorate a TrackingVolume
  virtual void
  decorate(TrackingVolume& volume) const = 0;
};

}  // namespace