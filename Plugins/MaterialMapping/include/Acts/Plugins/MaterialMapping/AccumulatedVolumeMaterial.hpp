// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedVolumeMaterial.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/MaterialMapping/AccumulatedMaterialProperties.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

// TODO: probably combine both classes into a single one
struct AccumulatedMaterial
{
   void
  accumulate(const Material& mat)
  {
	  m_eventX0 += mat.X0();
	  m_eventL0 += mat.L0();
	  m_eventA += mat.A();
	  m_eventZ += mat.Z();
	  m_eventRho += mat.rho();
	  m_totalEntries++;
	}

  std::pair<Material, unsigned int>
  totalAverage()
  {
	    if (m_totalEntires > 0) {
		double eventScalor = 1. / (double) m_totalEntries;
		m_totalX0 *= eventScalor;
		m_totalL0 *= eventScalor;
		m_totalA *= eventScalor;
		m_totalZ *= eventScalor;
		m_totalRho *= eventScalor;
		// Create the material
		Material mat(m_totalX0, m_totalL0, m_totalA, m_totalZ, m_totalRho);
		return std::pair<Material, unsigned int>(std::move(mat), m_totalEntries);
	}
	else{
	return std::pair<Material, unsigned int>(Material(), m_totalEntries);}
  }

private:
  float m_eventX0{0.};  //!< event: accumulate the contribution to X0
  float m_eventL0{0.};  //!< event: accumulate the contribution to L0
  float m_eventA{0.};         //!< event: accumulate the contribution to A
  float m_eventZ{0.};         //!< event: accumulate the contribution to Z
  float m_eventRho{0.};       //!< event: accumulate the contribution to rho

  unsigned int m_totalEntries{0};  //!< the number of events
};

/// @class AccumulatedSurfaceMaterial
///
/// This class is used by the SurfaceMaterialMapper in order to
/// accumulate/collect material information during the mapping process.
///
/// It performs event- and run-average when called, and returns
/// a new SurfaceMaterial object as a unique_ptr after finalisation
class AccumulatedVolumeMaterial
{
public:
  using AccumulatedVector = std::vector<AccumulatedMaterialProperties>;
  using AccumulatedMatrix = std::vector<AccumulatedVector>;

  /// Default Constructor - for homogeneous material
  ///
  /// @param splitFactor is the pre/post splitting directive
  AccumulatedSurfaceMaterial(double splitFactor = 0.);

  /// Explicit constructor with only full MaterialProperties,
  /// for one-dimensional binning.
  ///
  /// The split factors:
  ///    - 1. : oppositePre
  ///    - 0. : alongPre
  ///  ===> 1 Dimensional array
  ///
  /// @param binUtility defines the binning structure on the surface
  /// @param splitFactor is the pre/post splitting directive
  AccumulatedSurfaceMaterial(const BinUtility& binUtility,
                             double            splitFactor = 0.);

  /// Copy Constructor
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial(const AccumulatedSurfaceMaterial& asma) = default;

  /// Copy Move Constructor
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial(AccumulatedSurfaceMaterial&& asma) = default;

  /// Assignment Move operator
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial&
  operator=(AccumulatedSurfaceMaterial&& asma)
      = default;

  /// Assignment operator
  ///
  /// @param asma is the source object to be copied
  AccumulatedSurfaceMaterial&
  operator=(const AccumulatedSurfaceMaterial& asma)
      = default;

  /// Destructor
  ~AccumulatedSurfaceMaterial() = default;

  /// Return the BinUtility
  const BinUtility&
  binUtility() const;

  /// Assign a material properites object
  ///
  /// @param lp local position for the bin assignment
  /// @param mp material properties to be assigned
  void
  accumulate(const Vector2D&           lp,
             const MaterialProperties& mp,
             double                    pathCorrection = 1.);

  /// Assign a material properites object
  ///
  /// @param gp local position for the bin assignment
  /// @param mp material properties to be assigned
  void
  accumulate(const Vector3D&           gp,
             const MaterialProperties& mp,
             double                    pathCorrection = 1.);

  /// Assign a material properites object
  ///
  /// @param gp local position for the bin assignment
  /// @param mps Vector of recorded material properties to be assigned
  /// @param pathCorrection The geometric path correction due to incident
  void
  accumulate(const Vector3D& gp,
             const std::vector<std::pair<MaterialProperties, Vector3D>>& mps,
             double pathCorrection = 1.);

  /// Average the information accumulated during one event
  /// using the event weights
  void
  eventAverage();

  /// Total average creates SurfaceMaterial
  std::unique_ptr<const SurfaceMaterial>
  totalAverage();

  /// Access to the accumulated material
  const AccumulatedMatrix&
  accumulatedMaterial() const;

  /// Access to the split factor
  double
  splitFactor() const;

private:
  /// The helper for the bin finding
  BinUtility m_binUtility{};

  /// the split factor
  double m_splitFactor{0.};

  /// The stored accumulated material matrix
  AccumulatedMatrix m_accumulatedMaterial;
};

inline const BinUtility&
AccumulatedSurfaceMaterial::binUtility() const
{
  return (m_binUtility);
}

inline const AccumulatedSurfaceMaterial::AccumulatedMatrix&
AccumulatedSurfaceMaterial::accumulatedMaterial() const
{
  return (m_accumulatedMaterial);
}

inline double
AccumulatedSurfaceMaterial::splitFactor() const
{
  return m_splitFactor;
}
}