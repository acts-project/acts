// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <array>
#include <iosfwd>

namespace Acts {

/// Proxy marking a volume for material mapping, carrying three-axis binning.
///
/// Like ProtoSurfaceMaterial, the proxy carries an AxisSpec for each material
/// coordinate. Deferred ranges and boundary types are supplied by the mapping
/// consumer through MultiAxisSpec3D::buildMultiAxis. One bin in each direction
/// describes homogeneous material. The proxy always returns vacuum.
class ProtoVolumeMaterial final : public IVolumeMaterial {
 public:
  /// Construct a homogeneous proxy with three deferred single-bin axes.
  ProtoVolumeMaterial() = default;

  /// Construct from fully specified or deferred material binning.
  /// @param binning three-axis description of the future material map
  explicit ProtoVolumeMaterial(const MultiAxisSpec3D& binning);

  /// @return the material map's binning specification
  const MultiAxisSpec3D& binning() const { return m_binning; }

  /// Return an unbounded vacuum cell.
  /// @return the vacuum cell
  VolumeMaterialCell materialCell(const Vector3& /*position*/) const final {
    return VolumeMaterialCell(Material::Vacuum());
  }

  /// Write the proxy binning description.
  /// @param stream destination stream
  /// @return the modified stream
  std::ostream& toStream(std::ostream& stream) const final;

 private:
  MultiAxisSpec3D m_binning{std::array{AxisSpec::DeferredEquidistant(1),
                                       AxisSpec::DeferredEquidistant(1),
                                       AxisSpec::DeferredEquidistant(1)}};
};

}  // namespace Acts
