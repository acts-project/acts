// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

/// @class ProtoVolumeMaterial
///
/// @brief proxy to VolumeMaterial hand over BinUtility
///
/// The ProtoVolumeMaterial class acts as a proxy to the VolumeMaterial
/// to mark the volume on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity of
/// of the material map with the bin Utility.

class ProtoVolumeMaterial : public IVolumeMaterial {
 public:
  /// Constructor without BinUtility - homogenous material
  ProtoVolumeMaterial() = default;

  /// Constructor with BinUtility - multidimensional material
  ///
  /// @param binUtility a BinUtility determining the granularity
  ///        and binning of the material on the volume
  ProtoVolumeMaterial(const BinUtility& binUtility)
  : IVolumeMaterial(), m_binUtility(binUtility) {};

  /// Copy constuctor
  ///
  /// @param vmproxy The source proxy
  ProtoVolumeMaterial(const ProtoVolumeMaterial& vmproxy) = default;

  /// Copy move constuctor
  ///
  /// @param vmproxy The source proxy
  ProtoVolumeMaterial(ProtoVolumeMaterial&& vmproxy) = default;

  /// Destructor
  ///
  ~ProtoVolumeMaterial() override = default;

  /// Return the BinUtility
  const BinUtility& binUtility() const{
    return m_binUtility;
  }

  /// Assignment operator
  ///
  /// @param vmproxy The source proxy
  ProtoVolumeMaterial& operator=(const ProtoVolumeMaterial& vmproxy) = default;

  /// Return the material
  const Material& material(const Vector3D& /*position*/) const final {
    return m_material;
  }

 private:
  BinUtility m_binUtility;
  Material m_material;
};
}  // namespace Acts
