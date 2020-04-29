// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>

#include "Acts/Geometry/GeometryID.hpp"

namespace Acts {

class ISurfaceMaterial;
class IVolumeMaterial;

using SurfaceMaterialMap =
    std::map<GeometryID, std::shared_ptr<const ISurfaceMaterial>>;

using VolumeMaterialMap =
    std::map<GeometryID, std::shared_ptr<const IVolumeMaterial>>;

using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;
}  // namespace Acts

namespace FW {

/// @class IMaterialWriter
///
/// Interface definition for material writing
class IMaterialWriter {
 public:
  /// Virtual Destructor
  virtual ~IMaterialWriter() = default;

  /// The single wirter class
  ///
  /// @param detMaterial the detector material maps
  virtual void writeMaterial(const Acts::DetectorMaterialMaps& detMaterial) = 0;
};

/// @class MaterialWriterT
///
/// @tparam writer_t is the actual implementation
template <typename writer_t>
class MaterialWriterT : virtual public IMaterialWriter {
 public:
  /// Constructor
  ///
  /// @tparam writer_t the templated writer implementation
  ///
  /// @param impl the actaul implementation of the writer
  MaterialWriterT(writer_t impl) : m_impl(std::move(impl)) {}

  /// The single wirter class
  ///
  /// @param detMaterial the detector material maps
  void writeMaterial(const Acts::DetectorMaterialMaps& detMaterial) {
    m_impl.write(detMaterial);
  }

 private:
  /// The writer implementation
  writer_t m_impl;
};
}  // namespace FW
