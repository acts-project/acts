// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialRecord.h, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/MaterialMapping/AssignedMaterialSteps.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialStep.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

class Surface;

using RecordBin      = std::pair<MaterialProperties, size_t>;
using RecordVector   = std::vector<RecordBin>;
using MaterialRecord = std::vector<RecordVector>;

/// @class SurfaceMaterialRecord
///
/// @brief Records the material per layer during material mapping.
///
/// The Acts::SurfaceMaterialRecord class is used as a cache during the material
/// mapping process for the layers and hands back the final layer material.
/// It stores the accumulated material of a certain
/// layer in a matrix binned in a given Acts::BinUtility. Furthermore it also
/// stores a collection of all added material steps per track which can be used
/// to write out material maps of single layers.
///
/// The Acts::MaterialMapping class uses this class to add
/// material at a certain position on the layer which is transformed into the
/// corresponding bin of the grid. Furthermore it also uses it  to average the
/// material for each bin of the layer during the mapping process whenever
/// wanted
/// (e.g. after each run, after every event). In the end before handing
/// back the complete layer material an averaging needs to be done.
///

class SurfaceMaterialRecord
{
public:
  /// Default constructor - deleted
  SurfaceMaterialRecord() {}

  /// Constructor with BinUtility input
  /// @param surface is the according surface of this recrd
  /// @param binUtility for the the grid in which the material is binned on the
  /// layer
  SurfaceMaterialRecord(const Surface& surface, const BinUtility& binUtility);

  /// Default destructor
  ~SurfaceMaterialRecord() = default;

  /// Copy Constructor
  SurfaceMaterialRecord(const SurfaceMaterialRecord& lmrecord);

  /// Assignment operator
  SurfaceMaterialRecord&
  operator=(const SurfaceMaterialRecord& lmrecord);

  /// Adds MaterialProperties and weighs them over the steplength at a given
  /// position
  ///
  /// @param mStep material step
  /// @param pathCorrection is the path correction
  /// added at the given position
  void
  assignMaterialStep(const MaterialStep& mStep, double pathCorrection = 1.);

  /// Associate an empty step
  /// this is still needed, because the extrapolation
  /// might hit a layer, but no material to access was there
  ///
  /// @param mPosition is where the extrapolation
  ///    did hit this surface
  void
  assignEmptyStep(const Vector3D& mPosition);

  /// @return the surface pointer
  const Surface&
  surface() const;

  /// @return the bin utility
  const BinUtility&
  binUtility() const;

  /// this is the material without being averaged
  /// the averaging still has to be done in the mapper
  ///
  /// @return the full MaterialRecord
  const MaterialRecord&
  mappedMaterial() const;

private:
  /// remember the Surface
  const Surface* m_surface;

  /// two dimensional grid on which the material is binned
  std::unique_ptr<const BinUtility> m_binUtility;

  /// the material record
  MaterialRecord m_mappedMaterial;
};

inline const Surface&
SurfaceMaterialRecord::surface() const
{
  return (*m_surface);
}

inline const BinUtility&
SurfaceMaterialRecord::binUtility() const
{
  return (*m_binUtility);
}

inline const MaterialRecord&
SurfaceMaterialRecord::mappedMaterial() const
{
  return m_mappedMaterial;
}
}