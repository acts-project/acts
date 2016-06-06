// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PrismVolumeBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_VOLUMES_PRISMVOLUMEBOUNDS_H
#define ACTS_VOLUMES_PRISMVOLUMEBOUNDS_H 1

#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"

namespace Acts {

class Surface;
class PlaneSurface;
class TriangleBounds;

/**
 @class PrismVolumeBounds

 Bounds for the transcript of triangular prism

  BoundarySurfaceFace [index]:

      - negativeFaceXY     [0] : Triangular Acts::PlaneSurface,
                                 parallel to \f$ xy \f$ plane at negative \f$ z
 \f$
      - positiveFaceXY     [1] : Triangular Acts::PlaneSurface,
                                 parallel to \f$ xy \f$ plane at positive \f$ z
 \f$
      - face [2... n+1] : Rectangular  Acts::PlaneSurface

  */

class PrismVolumeBounds : public VolumeBounds
{
public:
  /**Default Constructor*/
  PrismVolumeBounds();

  /**Constructor - generic case (from float)*/
  PrismVolumeBounds(std::vector<std::pair<float, float>> xyvtx, float hlengthz);

  /**Constructor - generic case from (double)*/
  PrismVolumeBounds(std::vector<std::pair<double, double>> xyvtx,
                    double hlengthz);

  /**Copy Constructor */
  PrismVolumeBounds(const PrismVolumeBounds& bobo);

  /**Destructor */
  virtual ~PrismVolumeBounds();

  /**Assignment operator*/
  PrismVolumeBounds&
  operator=(const PrismVolumeBounds& bobo);

  /**Virtual constructor */
  PrismVolumeBounds*
  clone() const override;

  /**This method checks if position in the 3D volume frame is inside the
   * volume*/
  bool
  inside(const Vector3D&, double tol = 0.) const override;

  /** Method to decompose the Bounds into Surfaces */
  const std::vector<const Surface*>*
  decomposeToSurfaces(std::shared_ptr<Transform3D> transformPtr) const override;

  /**This method returns the set of xy generating vertices*/
  const std::vector<std::pair<TDD_real_t, TDD_real_t>>
  xyVertices() const;

  /**This method returns the halflength in local z*/
  double
  halflengthZ() const;

  /** Output Method for std::ostream */
  std::ostream&
  dump(std::ostream& sl) const override;

private:
  /** templated dump method */
  template <class T>
  T&
  dumpT(T& dt) const;

  /** method to construct side boundary planes */
  Acts::PlaneSurface*
  sideSurf(Transform3D, unsigned int, unsigned int) const;

  /** mirror the input vertices for down-side boundary */
  std::vector<std::pair<double, double>>
  mirror_xyVtx() const;

  /** assess ordering of vertices */
  int
  ordering() const;

  mutable std::vector<std::pair<TDD_real_t, TDD_real_t>>
             m_xyVtx;  //!< generating xy vertices
  TDD_real_t m_halfZ;  //!< halflength in z

  mutable Acts::TriangleBounds* m_baseBounds;  //!< base xy bounds
  mutable int                   m_ordering;    //!< cache vertex ordering
};

inline PrismVolumeBounds*
PrismVolumeBounds::clone() const
{
  return new PrismVolumeBounds(*this);
}

inline const std::vector<std::pair<TDD_real_t, TDD_real_t>>
PrismVolumeBounds::xyVertices() const
{
  return m_xyVtx;
}

inline double
PrismVolumeBounds::halflengthZ() const
{
  return m_halfZ;
}

template <class T>
T&
PrismVolumeBounds::dumpT(T& dt) const
{
  dt << std::setiosflags(std::ios::fixed);
  dt << std::setprecision(7);
  dt << "Acts::PrismVolumeBounds: (halfZ, generating vtx) = ";
  dt << "( " << m_halfZ << ")";
  for (unsigned int i = 0; i < m_xyVtx.size(); i++)
    dt << "(" << m_xyVtx.at(i).first << "," << m_xyVtx.at(i).second << ")";
  // return the modified stream
  return dt;
}
}

#endif  // ACTS_VOLUMES_PRISMVOLUMEBOUNDS_H
