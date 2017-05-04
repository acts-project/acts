// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialRecord.h, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIALPLUGINS_SURFACEMATERIALRECORD_H
#define ACTS_MATERIALPLUGINS_SURFACEMATERIALRECORD_H 1

#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Plugins/MaterialPlugins/AssignedMaterialSteps.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

namespace Acts {

class Surface;

typedef std::pair<MaterialProperties, size_t> RecordBin;
typedef std::vector<RecordBin>    RecordVector;
typedef std::vector<RecordVector> MaterialRecord;

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
  /// @param pos global position at which the material should be added
  /// @param materialSteps the material steps of a track which should be
  /// added at the given position
  void
  assignMaterialSteps(const AssignedMaterialSteps& aSteps);

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
  std::unique_ptr<BinUtility> m_binUtility;

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

#endif  // ACTS_MATERIALPLUGINS_SURFACEMATERIALRECORD_H
