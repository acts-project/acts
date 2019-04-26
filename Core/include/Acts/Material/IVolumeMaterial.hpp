// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IVolumeMaterial.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

namespace Acts {

class Material;

/// @class IVolumeMaterial
///
/// Virtual base class of volume material description
//
/// Material associated with a Volume (homogenous, binned, interpolated)
class IVolumeMaterial {
 public:
  /// Virtual Destructor
  virtual ~IVolumeMaterial() = default;

  /// Access to actual material
  ///
  /// @param position is the request position for the material call
  /// @todo interface to change including 'cell'
  virtual const Material& material(const Vector3D& position) const = 0;
};

}  // namespace Acts
