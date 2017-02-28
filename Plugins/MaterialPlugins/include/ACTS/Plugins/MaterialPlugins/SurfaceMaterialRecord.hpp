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

#ifndef ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H
#define ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H

#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

namespace Acts {

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
  /// Default constructor
  SurfaceMaterialRecord();
  
  /// Constructor with BinUtility input
  /// @param binutility the 2D grid in which the material is binned on the layer
  SurfaceMaterialRecord(const BinUtility* binutility);
  
  /// Default destructor
  ~SurfaceMaterialRecord() = default;
  
  /// Copy Constructor
  SurfaceMaterialRecord(const SurfaceMaterialRecord& lmrecord);
  
  /// Implicit contructor
  /// - uses the copy constructor
  SurfaceMaterialRecord*
  clone() const;
  
  /// Assignment operator
  SurfaceMaterialRecord&
  operator=(const SurfaceMaterialRecord& lmrecord);
  
  /// Adds MaterialProperties and weighs them over the steplength at a given
  /// position
  /// @param pos global position at which the material should be added
  /// @param layerMaterialSteps the material steps of a track which should be
  /// added at the given position
  void
  assignMaterialSteps(
      const Vector3D&                 pos,
      const std::vector<MaterialStep> materialSteps);
      
  /// @return returns all material steps per Track (original position) with
  /// their assigned material positions on the layer
  ///  @note this function is intended to be used to create material maps for
  ///  single layers if needed
  const std::vector<std::pair<const std::vector<MaterialStep>,
                              const Vector3D>>
  assignedMaterialSteps() const;
                              
  /// Possibility to average over the material given so far
  /// resets the sums and the counter of how often a certain bin was hit
  void
  averageMaterial();

  /// @return method for the final layer material
  /// given as a binned surface material
  std::shared_ptr<const Acts::BinnedSurfaceMaterial>
  surfaceMaterial() const;

private:
  /// two dimensional grid on which the material is binned
  const BinUtility* m_binUtility;
  
  /// two dimensional material matrix describing the material binned according
  /// to the binUtility
  std::vector<std::vector<const MaterialProperties*>> m_materialMatrix;
  
  /// the collection of assigned material steps per track, this collection
  std::vector<std::pair<const vector<MaterialStep>, const Vector3D>>
      m_matStepsAndAssignedPos;
};

inline const std::vector<std::pair<const std::vector<MaterialStep>,
                                   const Acts::Vector3D>>
SurfaceMaterialRecord::layerMaterialSteps() const
{
  return m_matStepsAndAssignedPos;
}
}

#endif  // ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H
