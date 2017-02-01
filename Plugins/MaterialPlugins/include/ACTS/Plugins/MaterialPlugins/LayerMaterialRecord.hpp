// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerMaterialRecord.h, ACTS project MaterialPlugins
///////////////////////////////////////////////////////////////////

#ifndef ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H
#define ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H

#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Utilities/BinUtility.hpp"

namespace Acts {

/// @class LayerMaterialRecord
///
/// @brief Records the material per layer during material mapping.
///
/// The Acts::LayerMaterialRecord class is used as a cache during the material
/// mapping process for the layers and hands back the final layer material.
/// The user will never directly inderact with this class.
/// It stores the current material of a certain
/// layer in a matrix binned in a given Acts::BinUtility. It is possible to add
/// material at a certain position on the layer which is transformed into the
/// corresponding bin of the grid. It is possible to average the material for
/// each bin of the layer during the mapping process whenever wanted
/// (e.g. after each run, after every event). In the end before handing
/// back the complete layer material an averaging needs to be done.
///

class LayerMaterialRecord
{
public:
  /// Default constructor
  LayerMaterialRecord();
  /// Constructor with BinUtility input
  /// @param binutility the 2D grid in which the material is binned on the layer
  LayerMaterialRecord(const BinUtility* binutility);
  /// Default destructor
  ~LayerMaterialRecord() = default;
  /// Copy Constructor
  LayerMaterialRecord(const LayerMaterialRecord& lmrecord);
  /// Implicit contructor
  /// - uses the copy constructor
  LayerMaterialRecord*
  clone() const;
  /// Assignment operator
  LayerMaterialRecord&
  operator=(const LayerMaterialRecord& lmrecord);
  /// Adds MaterialProperties and weighs them over the steplength at a given
  /// position
  /// @param pos global position at which the material should be added
  /// @param newMaterial the material properties (material + step length) at the
  /// given position
  void
  addLayerMaterialProperties(const Acts::Vector3D&           pos,
                             const Acts::MaterialProperties* newMaterial);
  /// Possibility to average over the material given so far
  /// resets the sums and the counter of how often a certain bin was hit
  void
  averageMaterial();
  /// @return method for the final layer material
  /// given as a binned surface material
  std::shared_ptr<const Acts::BinnedSurfaceMaterial>
  layerMaterial() const;

private:
  /// two dimensional grid on which the material is binned
  const BinUtility* m_binUtility;
  /// two dimensional material matrix describing the material binned according
  /// to the binUtility
  std::vector<std::vector<const Acts::MaterialProperties*>> m_materialMatrix;
};
}

#endif  // ACTS_MATERIALPLUGINS_LAYERMATERIALRECORD_H
