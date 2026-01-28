// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <sstream>

namespace Acts {

/// This enum describes the type of surface material mapping
enum class MappingType : std::int8_t {
  PreMapping = -1,
  Default = 0,
  PostMapping = 1,
  Sensor = 2
};

/// @ingroup material
///
/// Base class of all surface-based material description
///
/// The class supplies references to @ref MaterialSlab that are associated to a
/// surface, extended by certain special representations (binned, homogeneous).
/// The concrete @ref MaterialSlab can depend on the local position on the
/// surface.
///
class ISurfaceMaterial {
 public:
  /// Constructor
  ISurfaceMaterial() = default;

  /// Constructor
  ///
  /// @param splitFactor is the splitting ratio between pre/post update
  explicit ISurfaceMaterial(double splitFactor) : m_splitFactor(splitFactor) {}

  /// Constructor
  ///
  /// @param splitFactor is the splitting ratio between pre/post update
  /// @param mappingType is the type of surface mapping associated to the surface
  explicit ISurfaceMaterial(double splitFactor, MappingType mappingType)
      : m_splitFactor(splitFactor), m_mappingType(mappingType) {}

  /// Destructor
  virtual ~ISurfaceMaterial() = default;

  /// Scale material
  ///
  /// @param factor is the scale factor applied
  /// @return Reference to this material object for chaining
  virtual ISurfaceMaterial& scale(double factor) = 0;

  /// Return method for full material description of the Surface
  /// - from local coordinate on the surface
  ///
  /// @param lp is the local position used for the (eventual) lookup
  ///
  /// @return const MaterialSlab
  virtual const MaterialSlab& materialSlab(const Vector2& lp) const = 0;

  /// Return method for full material description of the Surface
  /// - from the global coordinates
  ///
  /// @param gp is the global position used for the (eventual) lookup
  ///
  /// @return const MaterialSlab
  virtual const MaterialSlab& materialSlab(const Vector3& gp) const = 0;

  /// Update pre factor
  ///
  /// @param pDir is the positive direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  /// @return Factor for material scaling based on direction and update stage
  double factor(Direction pDir, MaterialUpdateStage mStage) const;

  /// Return the type of surface material mapping
  ///
  /// @return The mapping type indicating how material is associated with the surface
  MappingType mappingType() const { return m_mappingType; }

  /// Return method for fully scaled material description of the Surface
  /// - from local coordinate on the surface
  ///
  /// @param lp is the local position used for the (eventual) lookup
  /// @param pDir is the positive direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  ///
  /// @return MaterialSlab
  virtual MaterialSlab materialSlab(const Vector2& lp, Direction pDir,
                                    MaterialUpdateStage mStage) const;

  /// Return method for full material description of the Surface
  /// - from the global coordinates
  ///
  /// @param gp is the global position used for the (eventual) lookup
  /// @param pDir is the positive direction through the surface
  /// @param mStage is the material update directive (onapproach, full, onleave)
  ///
  /// @return MaterialSlab
  virtual MaterialSlab materialSlab(const Vector3& gp, Direction pDir,
                                    MaterialUpdateStage mStage) const;

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual ISurfaceMaterial::toStream method
  ///
  /// @return modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const ISurfaceMaterial& sm) {
    sm.toStream(out);
    return out;
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param sl Output stream to write to
  /// @return Reference to the output stream for chaining
  virtual std::ostream& toStream(std::ostream& sl) const = 0;

  /// @brief output into a string
  ///
  /// @return the string representation
  std::string toString() const {
    std::stringstream sstrm;
    toStream(sstrm);
    return sstrm.str();
  }

 protected:
  /// the split factor in favour of oppositePre
  double m_splitFactor{1.};

  /// Use the default mapping type by default
  MappingType m_mappingType{MappingType::Default};
};

inline double ISurfaceMaterial::factor(Direction pDir,
                                       MaterialUpdateStage mStage) const {
  if (mStage == Acts::MaterialUpdateStage::FullUpdate) {
    return 1.;
  } else if (mStage == Acts::MaterialUpdateStage::PreUpdate) {
    return pDir == Direction::Negative() ? m_splitFactor : 1 - m_splitFactor;
  } else /*if (mStage == Acts::MaterialUpdateStage::PostUpdate)*/ {
    return pDir == Direction::Positive() ? m_splitFactor : 1 - m_splitFactor;
  }
}

inline MaterialSlab ISurfaceMaterial::materialSlab(
    const Vector2& lp, Direction pDir, MaterialUpdateStage mStage) const {
  // The plain material properties associated to this bin
  MaterialSlab plainMatProp = materialSlab(lp);
  // Scale if you have material to scale
  if (!plainMatProp.isVacuum()) {
    double scaleFactor = factor(pDir, mStage);
    if (scaleFactor == 0.) {
      return MaterialSlab::Nothing();
    }
    plainMatProp.scaleThickness(static_cast<float>(scaleFactor));
  }
  return plainMatProp;
}

inline MaterialSlab ISurfaceMaterial::materialSlab(
    const Vector3& gp, Direction pDir, MaterialUpdateStage mStage) const {
  // The plain material properties associated to this bin
  MaterialSlab plainMatProp = materialSlab(gp);
  // Scale if you have material to scale
  if (!plainMatProp.isVacuum()) {
    double scaleFactor = factor(pDir, mStage);
    if (scaleFactor == 0.) {
      return MaterialSlab::Nothing();
    }
    plainMatProp.scaleThickness(static_cast<float>(scaleFactor));
  }
  return plainMatProp;
}

}  // namespace Acts
