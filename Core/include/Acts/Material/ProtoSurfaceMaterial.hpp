// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <iosfwd>
#include <utility>
#include <vector>

namespace Acts {

/// @addtogroup material
/// @{

///
/// @brief proxy to SurfaceMaterial carrying directed proto axis binning
/// and a global-to-local transform
///
/// The ProtoSurfaceMaterial class acts as a proxy to the SurfaceMaterial
/// to mark the layers and surfaces on which the material should be mapped on
/// at construction time of the geometry and to hand over the granularity of
/// of the material map with directed proto axes.
class ProtoSurfaceMaterial : public ISurfaceMaterial {
 public:
  /// Constructor without binningType - homogeneous material
  ProtoSurfaceMaterial() = default;

  /// Constructor with directed proto axes and global-to-local transform
  /// @param directedProtoAxes axis description for the material map binning
  /// @param globalToLocalTransform transform from global to local 3D frame
  /// @param mappingType is the type of surface mapping associated to the surface
  explicit ProtoSurfaceMaterial(
      std::vector<DirectedProtoAxis> directedProtoAxes,
      Transform3 globalToLocalTransform = Transform3::Identity(),
      MappingType mappingType = MappingType::Default)
      : ISurfaceMaterial(1., mappingType),
        m_directedProtoAxes(std::move(directedProtoAxes)),
        m_globalToLocalTransform(std::move(globalToLocalTransform)) {}

  /// Copy constructor
  ///
  /// @param smproxy The source proxy
  ProtoSurfaceMaterial(const ProtoSurfaceMaterial& smproxy) = default;

  /// Copy move constructor
  ///
  /// @param smproxy The source proxy
  ProtoSurfaceMaterial(ProtoSurfaceMaterial&& smproxy) noexcept = default;

  /// Destructor
  ~ProtoSurfaceMaterial() override = default;

  /// Assignment operator
  ///
  /// @param smproxy The source proxy
  /// @return Reference to this object
  ProtoSurfaceMaterial& operator=(const ProtoSurfaceMaterial& smproxy) =
      default;

  /// Assignment move operator
  ///
  /// @param smproxy The source proxy
  /// @return Reference to this object
  ProtoSurfaceMaterial& operator=(ProtoSurfaceMaterial&& smproxy) noexcept =
      default;

  /// Scale operation - dummy implementation
  ///
  /// @return Reference to this object
  ProtoSurfaceMaterial& scale(double /*factor*/) final { return (*this); }

  /// Return the directed proto axes
  /// @return Reference to the binning axis descriptors
  const std::vector<DirectedProtoAxis>& directedProtoAxes() const {
    return m_directedProtoAxes;
  }

  /// Return transform from global to local 3D frame
  /// @return Reference to transform
  const Transform3& globalToLocalTransform() const {
    return m_globalToLocalTransform;
  }

  /// Return a BinUtility representation of this proxy (on-the-fly)
  [[deprecated(
      "ProtoSurfaceMaterial::binning() is deprecated. "
      "Use directedProtoAxes() and globalToLocalTransform() instead.")]]
  BinUtility binning() const {
    BinUtility converted(m_globalToLocalTransform.inverse());
    for (const auto& directedProtoAxis : m_directedProtoAxes) {
      converted += BinUtility(BinningData(directedProtoAxis));
    }
    return converted;
  }

  /// Return method for full material description of the Surface - from local
  /// coordinates
  ///
  /// @return will return dummy material
  const MaterialSlab& materialSlab(const Vector2& /*lp*/) const final {
    return (m_materialSlab);
  }

  /// Return method for full material description of the Surface - from the
  /// global coordinates
  ///
  /// @return will return dummy material
  const MaterialSlab& materialSlab(const Vector3& /*gp*/) const final {
    return (m_materialSlab);
  }

  using ISurfaceMaterial::materialSlab;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the output stream
  /// @return The output stream
  std::ostream& toStream(std::ostream& sl) const final {
    sl << "Acts::ProtoSurfaceMaterial : " << std::endl;
    sl << m_directedProtoAxes << std::endl;
    return sl;
  }

 private:
  /// Directed axis descriptions.
  std::vector<DirectedProtoAxis> m_directedProtoAxes;

  /// Transform from global to local 3D frame.
  Transform3 m_globalToLocalTransform = Transform3::Identity();

  /// Dummy material properties
  MaterialSlab m_materialSlab = MaterialSlab::Nothing();
};

/// @}

}  // namespace Acts
