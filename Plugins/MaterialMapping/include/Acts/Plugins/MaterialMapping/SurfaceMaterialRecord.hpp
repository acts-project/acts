// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceRecordMatrix.h, Acts project MaterialPlugins
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/MaterialMapping/AssignedMaterialSteps.hpp"
#include "Acts/Plugins/MaterialMapping/MaterialStep.hpp"
#include "Acts/Utilities/BinUtility.hpp"

namespace Acts {

class Surface;

using RecordBin    = std::pair<MaterialProperties, size_t>;
using RecordVector = std::vector<RecordBin>;
using RecordMatrix = std::vector<RecordVector>;

/// @class SurfaceRecordMatrix
///
/// @brief Records the material per layer during material mapping.
///
/// The Acts::SurfaceRecordMatrix class is used as a cache during the material
/// mapping process for the layers and hands back the final layer material.
/// It stores the accumulated material of a certain layer in a matrix binned
/// in a given Acts::BinUtility. Furthermore it also stores a collection of all
/// added material steps per track which can be used to write out material maps
/// of single layers.
///
/// The MaterialMapping class uses this class to add material at a certain
/// position on the layer which is transformed into the corresponding bin of
/// the grid. Furthermore it also uses it  to average the material for each bin
/// of the layer during the mapping process whenever wanted (e.g. after each
/// run, after every event). In the end before handing back the complete layer
/// material an averaging needs to be done.

class SurfaceRecordMatrix
{
public:
  /// Default constructor - default
  SurfaceRecordMatrix() = default;

  /// Constructor with BinUtility input
  /// @param surface The according surface of this record
  /// @param binUtility describes the binning on the surface
  SurfaceRecordMatrix(const Surface& surface, const BinUtility& binUtility);

  /// Default destructor
  ~SurfaceRecordMatrix() = default;

  /// Copy Constructor
  ///
  /// @param lmrecord The source matrix record
  SurfaceRecordMatrix(const SurfaceRecordMatrix& lmrecord);

  /// Assignment operator
  ///
  /// @param lmrecord The source matrix record
  SurfaceRecordMatrix&
  operator=(const SurfaceRecordMatrix& lmrecord);

  /// Adds MaterialProperties and weights them over the steplength
  /// at a given position
  ///
  /// @param mStep material step
  /// @param pathCorrection The incidence angle correction
  void
  assignMaterialStep(const MaterialStep& mStep, double pathCorrection = 1.);

  /// Associate an empty step this is still needed, because the extrapolation
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
  /// @return the full RecordMatrix
  const RecordMatrix&
  mappedMaterial() const;

private:
  /// remember the Surface
  const Surface* m_surface;

  /// two dimensional grid on which the material is binned
  std::unique_ptr<const BinUtility> m_binUtility;

  /// the material record
  RecordMatrix m_mappedMaterial;
};

inline const Surface&
SurfaceRecordMatrix::surface() const
{
  return (*m_surface);
}

inline const BinUtility&
SurfaceRecordMatrix::binUtility() const
{
  return (*m_binUtility);
}

inline const RecordMatrix&
SurfaceRecordMatrix::mappedMaterial() const
{
  return m_mappedMaterial;
}
}