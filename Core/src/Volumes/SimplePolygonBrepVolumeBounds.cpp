// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SimplePolygonBrepVolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/SimplePolygonBrepVolumeBounds.hpp"
#include "ACTS/Surfaces/PlaneSurface.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Surfaces/SubtractedPlaneSurface.hpp"
#include "ACTS/Volumes/CombinedVolumeBounds.hpp"
#include "ACTS/Volumes/CuboidVolumeBounds.hpp"
#include "ACTS/Volumes/PrismVolumeBounds.hpp"
#include "ACTS/Volumes/Volume.hpp"
#include "ACTS/Volumes/VolumeExcluder.hpp"
// STD/STL
#include <iomanip>
#include <iostream>
#include <math.h>

Acts::SimplePolygonBrepVolumeBounds::SimplePolygonBrepVolumeBounds()
  : VolumeBounds()
  , m_halfX(0.)
  , m_halfY(0.)
  , m_halfZ(0.)
  , m_ordering(-1)
  , m_combinedVolume(0)
  , m_envelope(0)
{
  //@TODO an object created by the default constructor cannot be copied or
  //assigned.
}

Acts::SimplePolygonBrepVolumeBounds::SimplePolygonBrepVolumeBounds(
    std::vector<std::pair<float, float>> xyVtx,
    float halez)
  : VolumeBounds()
  , m_halfX(0.)
  , m_halfY(0.)
  , m_halfZ(halez)
  , m_ordering(-1)
  , m_combinedVolume(0)
  , m_envelope(0)
{
  m_xyVtx.resize(xyVtx.size());
  double xmin = xyVtx.at(0).first;
  double xmax = xyVtx.at(0).first;
  double ymin = xyVtx.at(0).second;
  double ymax = xyVtx.at(0).second;
  for (unsigned int i = 0; i < xyVtx.size(); i++) {
    m_xyVtx.at(i)                       = xyVtx.at(i);
    if (xyVtx.at(i).first < xmin) xmin  = xyVtx.at(i).first;
    if (xyVtx.at(i).first > xmax) xmax  = xyVtx.at(i).first;
    if (xyVtx.at(i).second < ymin) ymin = xyVtx.at(i).second;
    if (xyVtx.at(i).second > ymax) ymax = xyVtx.at(i).second;
  }
  double ehalfX = 0.5 * (xmax - xmin);
  double ehalfY = 0.5 * (ymax - ymin);
  m_halfX       = fmax(fabs(xmax), fabs(xmin));
  m_halfY       = fmax(fabs(ymax), fabs(ymin));
  Acts::Transform3D transXY(
      Acts::Translation3D(Acts::Vector3D(0.5 * (xmin + xmax), 0., 0.))
      * Acts::Translation3D(Acts::Vector3D(0., 0.5 * (ymin + ymax), 0.)));
  m_envelope
      = new Acts::Volume(new Acts::Transform3D(transXY),
                         new Acts::CuboidVolumeBounds(ehalfX, ehalfY, m_halfZ));

  processSubVols();
}

Acts::SimplePolygonBrepVolumeBounds::SimplePolygonBrepVolumeBounds(
    std::vector<std::pair<double, double>> xyVtx,
    double halez)
  : VolumeBounds()
  , m_halfX(0.)
  , m_halfY(0.)
  , m_halfZ(halez)
  , m_ordering(-1)
  , m_combinedVolume(0)
  , m_envelope(0)
{
  m_xyVtx.resize(xyVtx.size());
  double xmin = xyVtx.at(0).first;
  double xmax = xyVtx.at(0).first;
  double ymin = xyVtx.at(0).second;
  double ymax = xyVtx.at(0).second;
  for (unsigned int i = 0; i < xyVtx.size(); i++) {
    m_xyVtx.at(i)                       = xyVtx.at(i);
    if (xyVtx.at(i).first < xmin) xmin  = xyVtx.at(i).first;
    if (xyVtx.at(i).first > xmax) xmax  = xyVtx.at(i).first;
    if (xyVtx.at(i).second < ymin) ymin = xyVtx.at(i).second;
    if (xyVtx.at(i).second > ymax) ymax = xyVtx.at(i).second;
  }
  double ehalfX = 0.5 * (xmax - xmin);
  double ehalfY = 0.5 * (ymax - ymin);
  m_halfX       = fmax(fabs(xmax), fabs(xmin));
  m_halfY       = fmax(fabs(ymax), fabs(ymin));
  Acts::Transform3D transXY(
      Acts::Translation3D(Acts::Vector3D(0.5 * (xmin + xmax), 0., 0.))
      * Acts::Translation3D(Acts::Vector3D(0., 0.5 * (ymin + ymax), 0.)));
  m_envelope
      = new Acts::Volume(new Acts::Transform3D(transXY),
                         new Acts::CuboidVolumeBounds(ehalfX, ehalfY, m_halfZ));

  processSubVols();
}

Acts::SimplePolygonBrepVolumeBounds::SimplePolygonBrepVolumeBounds(
    const Acts::SimplePolygonBrepVolumeBounds& trabo)
  : VolumeBounds()
  , m_halfX(trabo.m_halfX)
  , m_halfY(trabo.m_halfY)
  , m_halfZ(trabo.m_halfZ)
  , m_ordering(trabo.m_ordering)
  , m_combinedVolume(trabo.m_combinedVolume->clone())
  , m_envelope(trabo.m_envelope->clone())
{
  m_xyVtx.resize(trabo.m_xyVtx.size());
  for (unsigned int i = 0; i < m_xyVtx.size(); i++)
    m_xyVtx.at(i)     = trabo.m_xyVtx.at(i);
}

Acts::SimplePolygonBrepVolumeBounds::~SimplePolygonBrepVolumeBounds()
{
  delete m_combinedVolume;
  delete m_envelope;
}

Acts::SimplePolygonBrepVolumeBounds&
Acts::SimplePolygonBrepVolumeBounds::
operator=(const Acts::SimplePolygonBrepVolumeBounds& trabo)
{
  if (this != &trabo) {
    m_halfX = trabo.m_halfX;
    m_halfY = trabo.m_halfY;
    m_halfZ = trabo.m_halfZ;
    m_xyVtx.resize(trabo.m_xyVtx.size());
    for (unsigned int i = 0; i < m_xyVtx.size(); i++)
      m_xyVtx.at(i)     = trabo.m_xyVtx.at(i);
    delete m_combinedVolume;
    delete m_envelope;
    m_combinedVolume = trabo.m_combinedVolume->clone();
    m_envelope       = trabo.m_envelope->clone();
  }
  return *this;
}

const std::vector<const Acts::Surface*>*
Acts::SimplePolygonBrepVolumeBounds::decomposeToSurfaces(
    std::shared_ptr<Acts::Transform3D> transformPtr) const
{
  std::vector<const Acts::Surface*>* retsf
      = new std::vector<const Acts::Surface*>;
  // the transform
  Acts::Transform3D transform = (transformPtr == nullptr)
      ? Acts::Transform3D::Identity()
      : (*transformPtr.get());
  Acts::Transform3D* tTransform = 0;

  // face surfaces xy
  //  (1) - at negative local z
  tTransform = new Acts::Transform3D(
      transform * Acts::Translation3D(Acts::Vector3D(0., 0., -m_halfZ)));
  Acts::PlaneSurface xymPlane(std::shared_ptr<Acts::Transform3D>(tTransform),
                              new Acts::RectangleBounds(m_halfX, m_halfY));
  Acts::VolumeExcluder* volExcl = new Acts::VolumeExcluder(
      new Acts::Volume(*m_combinedVolume,
                       Acts::Transform3D(Acts::Translation3D(
                           Acts::Vector3D(0., 0., -m_halfZ)))));
  retsf->push_back(new Acts::SubtractedPlaneSurface(xymPlane, volExcl, true));
  //  (2) - at positive local z
  tTransform = new Acts::Transform3D(
      transform * Acts::Translation3D(Acts::Vector3D(0., 0., m_halfZ)));
  Acts::PlaneSurface xyPlane(std::shared_ptr<Acts::Transform3D>(tTransform),
                             new Acts::RectangleBounds(m_halfX, m_halfY));
  volExcl = new Acts::VolumeExcluder(new Acts::Volume(
      *m_combinedVolume,
      Acts::Transform3D(Acts::Translation3D(Acts::Vector3D(0., 0., m_halfZ)))));
  retsf->push_back(new Acts::SubtractedPlaneSurface(xyPlane, volExcl, true));
  // loop over xy vertices
  //  (3)
  for (unsigned int iv = 0; iv < m_xyVtx.size(); iv++) {
    if (iv != m_xyVtx.size() - 1)
      retsf->push_back(sideSurf(transform, iv, iv + 1));
    else
      retsf->push_back(sideSurf(transform, iv, 0));
  }

  return retsf;
}

// faces in xy
Acts::PlaneSurface*
Acts::SimplePolygonBrepVolumeBounds::sideSurf(Acts::Transform3D transform,
                                              unsigned int      iv1,
                                              unsigned int      iv2) const
{
  Acts::PlaneSurface* plane = 0;

  double xdif  = m_xyVtx.at(iv2).first - m_xyVtx.at(iv1).first;
  double ydif  = m_xyVtx.at(iv2).second - m_xyVtx.at(iv1).second;
  double xsize = sqrt(xdif * xdif + ydif * ydif);

  double ori = m_ordering > 0 ? -1. : 1.;

  Acts::Vector3D pos(0.5 * (m_xyVtx.at(iv1).first + m_xyVtx.at(iv2).first),
                     0.5 * (m_xyVtx.at(iv1).second + m_xyVtx.at(iv2).second),
                     0.);
  double phi                   = ori * ydif < 0 ? M_PI / 2 : -M_PI / 2;
  if (ori > 0 && ydif > 0) phi = M_PI / 2;
  if (fabs(xdif) > 1e-6) {
    phi = atan(ydif / xdif);
    if (xdif < 0) phi += M_PI;
  }

  Acts::Transform3D* tTransform = new Acts::Transform3D(
      transform * Acts::Translation3D(pos)
      * Acts::AngleAxis3D(phi, Acts::Vector3D(0., 0., 1.))
      * Acts::AngleAxis3D(-ori * 90 * M_PI / 180, Acts::Vector3D(1., 0., 0.)));
  plane
      = new Acts::PlaneSurface(std::shared_ptr<Acts::Transform3D>(tTransform),
                               new Acts::RectangleBounds(0.5 * xsize, m_halfZ));

  return plane;
}

bool
Acts::SimplePolygonBrepVolumeBounds::inside(const Acts::Vector3D& pos,
                                            double                tol) const
{
  return (m_combinedVolume->inside(pos, tol));
}

void
Acts::SimplePolygonBrepVolumeBounds::processSubVols() const
{
  // translate into prisms (triangulate)
  Acts::Volume* cVol = 0;
#ifdef TRKDETDESCR_USEFLOATPRECISON
#define double float
#endif
  std::vector<std::pair<double, double>> triangles = TriangulatePolygonCheck(
      m_xyVtx);  //@TODO change argument to const vector<pair< > >
  std::vector<std::pair<double, double>> vertices;
#ifdef TRKDETDESCR_USEFLOATPRECISON
#undef double
#endif
  for (unsigned int i = 0; i < triangles.size(); i = i + 3) {
    vertices.push_back(triangles.at(i));
    vertices.push_back(triangles.at(i + 1));
    vertices.push_back(triangles.at(i + 2));
    Acts::Volume* newVol
        = new Acts::Volume(0, new Acts::PrismVolumeBounds(vertices, m_halfZ));
    if (cVol)
      cVol = new Acts::Volume(
          0, new Acts::CombinedVolumeBounds(cVol, newVol, false));
    else
      cVol = newVol;
    vertices.clear();
  }
  m_combinedVolume = cVol;
}

// ostream operator overload
std::ostream&
Acts::SimplePolygonBrepVolumeBounds::dump(std::ostream& sl) const
{
  return dumpT<std::ostream>(sl);
}

//////////////////////////////////////////////////////////////////////////
// Triangulate Polygon
// M. Wolter
//////////////////////////////////////////////////////////////////////////

bool
Acts::SimplePolygonBrepVolumeBounds::Xor(bool x, bool y) const
// XOR: Arguments are negated to ensure that they are 0/1. Then the bitwise Xor
// operator may
// apply.
{
  return !x ^ !y;
}

bool
Acts::SimplePolygonBrepVolumeBounds::Left(
    std::pair<TDD_real_t, TDD_real_t> a,
    std::pair<TDD_real_t, TDD_real_t> b,
    std::pair<TDD_real_t, TDD_real_t> c) const
// Returns true iff c is strictly to the left of the directed line through a to
// b.
{
  double CrossZ = (b.first - a.first) * (c.second - a.second)
      - (c.first - a.first) * (b.second - a.second);
  if (m_ordering == 1) return (CrossZ >= 0.);
  if (m_ordering == 0) return (CrossZ < 0.);
  return false;
}

bool
Acts::SimplePolygonBrepVolumeBounds::Intersect(
    std::pair<TDD_real_t, TDD_real_t> a,
    std::pair<TDD_real_t, TDD_real_t> b,
    std::pair<TDD_real_t, TDD_real_t> c,
    std::pair<TDD_real_t, TDD_real_t> d) const
// Returns true iff segments ab and cd intersect
{
  return Xor(Left(a, b, c), Left(a, b, d)) && Xor(Left(c, d, a), Left(c, d, b));
}

bool
Acts::SimplePolygonBrepVolumeBounds::InCone(
    int i,
    int j,
    std::vector<std::pair<TDD_real_t, TDD_real_t>> inputVertices) const
// 	Returns true iff the diagonal (i,j) is internal to the polygon in
//  the neighborhood of the i endpoint.
{
  int iPlus1  = (i + 1) % inputVertices.size();
  int iMinus1 = (i + inputVertices.size() - 1) % inputVertices.size();

  /* If P[i] is a convex vertex [ i+1 left or on (i-1,i) ]. */
  if (Left(inputVertices.at(iMinus1),
           inputVertices.at(i),
           inputVertices.at(iPlus1)))
    return Left(inputVertices.at(i),
                inputVertices.at(j),
                inputVertices.at(iMinus1))
        && Left(inputVertices.at(j),
                inputVertices.at(i),
                inputVertices.at(iPlus1));

  /* Assume (i-1,i,i+1) not collinear. */
  /* else v_i is reflex. */
  else
    return !(
        Left(inputVertices.at(i), inputVertices.at(j), inputVertices.at(iPlus1))
        && Left(inputVertices.at(j),
                inputVertices.at(i),
                inputVertices.at(iMinus1)));
}

bool
Acts::SimplePolygonBrepVolumeBounds::Diagonalie(
    int i,
    int j,
    std::vector<std::pair<TDD_real_t, TDD_real_t>> inputVertices) const
{
  // Returns TRUE iff (v_i, v_j) is a proper internal *or* external diagonal of
  // this polygon,
  // *ignoring edges incident to v_i and v_j*.

  /* For each edge (k,k+1) of P */
  for (int k = 0; k < (int)inputVertices.size(); k++) {
    int kPlus1 = (k + 1) % inputVertices.size();

    /* Skip edges incident to i or j */
    if (!((k == i) || (kPlus1 == i) || (k == j) || (kPlus1 == j)))
      if (Intersect(inputVertices.at(i),
                    inputVertices.at(j),
                    inputVertices.at(k),
                    inputVertices.at(kPlus1)))
        return false;
  }
  return true;
}

bool
Acts::SimplePolygonBrepVolumeBounds::Diagonal(
    int i,
    int j,
    std::vector<std::pair<TDD_real_t, TDD_real_t>> inputVertices) const
// Returns TRUE iff (v_i, v_j) is a proper internal diagonal of P.
{
  // std::cout<<"MW Diagonal "<<i<<" "<<j<<" "<<InCone(i,j, inputVertices)<<"
  // "<<Diagonalie(i,j, inputVertices)<<std::endl;
  return InCone(i, j, inputVertices) && Diagonalie(i, j, inputVertices);
}

std::vector<std::pair<TDD_real_t, TDD_real_t>>
Acts::SimplePolygonBrepVolumeBounds::TriangulatePolygon(
    const std::vector<std::pair<TDD_real_t, TDD_real_t>>& Vertices) const
{
  // Subtracting ears method
  //
  // One way to triangulate a simple polygon is by using the assertion that any
  // simple polygon without holes
  //  has at least two so called 'ears'. An ear is a triangle with two sides on
  //  the edge of the polygon and the
  //  other one completely inside it. The algorithm then consists of finding
  //  such an ear, removing it from the
  //  polygon (which results in a new polygon that still meets the conditions)
  //  and repeating until there is
  //  only one triangle left.
  //

  int NSize = Vertices.size();
  std::vector<std::pair<TDD_real_t, TDD_real_t>> outTriangles;
  std::vector<std::pair<TDD_real_t, TDD_real_t>> inputVertices;
  for (int i = 0; i < NSize; i++) inputVertices.push_back((Vertices).at(i));

  // for (int i; i<NSize;i++) std::cout<<"MW input vertices:
  // "<<inputVertices.at(i).first<<" "<<inputVertices.at(i).second<<std::endl;
  // Triangulates this polygon and saves triangle edges in TriPoly.
  // Triangles are stored CCW, with each set of 3 consecutive points in TriPoly
  // representing 1 triangle.
  // Assumes this polygon is closed.

  if (NSize < 4) return inputVertices;

  // Start triangulating
  int VerticesLeft = NSize;
  while (VerticesLeft > 3) {
    // std::cout<<"MW vertices left "<<VerticesLeft<<std::endl;
    bool bCornerCut = false;
    for (int i = 0; i < VerticesLeft; i++) {
      int iPlus1                         = i + 1;
      if (iPlus1 == VerticesLeft) iPlus1 = 0;
      int iPlus2                         = (iPlus1 + 1);
      if (iPlus2 == VerticesLeft) iPlus2 = 0;

      if (Diagonal(i, iPlus2, inputVertices)) {
        outTriangles.push_back(inputVertices.at(i));
        outTriangles.push_back(inputVertices.at(iPlus1));
        outTriangles.push_back(inputVertices.at(iPlus2));

        inputVertices.erase(inputVertices.begin() + iPlus1);
        VerticesLeft--;
        bCornerCut = true;
        break;
      }
    }
    if (!bCornerCut) {  // Error - bad poly
      // std::cout<<"MW 	Error - bad poly"<<std::endl;
      std::vector<std::pair<TDD_real_t, TDD_real_t>> out;
      return out;
    }
  }

  if (VerticesLeft == 3) {
    outTriangles.push_back(inputVertices.at(0));
    outTriangles.push_back(inputVertices.at(1));
    outTriangles.push_back(inputVertices.at(2));
    inputVertices.erase(inputVertices.begin() + 1);
    VerticesLeft--;
  }

  return outTriangles;
}

std::vector<std::pair<TDD_real_t, TDD_real_t>>
Acts::SimplePolygonBrepVolumeBounds::TriangulatePolygonCheck(
    const std::vector<std::pair<TDD_real_t, TDD_real_t>>& Vertices) const
{
  // Perform triangulation. Check the orientation of the verices in the polygon
  // m_ordering   = -1    not set
  // m_ordering   =  1    anticlockwise
  // m_ordering   =  0    clockwise

  if (m_ordering == -1) m_ordering = 1;
  std::vector<std::pair<double, double>> outTriangles
      = TriangulatePolygon(Vertices);
  if (outTriangles.size() == 0) {
    m_ordering   = -m_ordering + 1;
    outTriangles = TriangulatePolygon(Vertices);
  }

  return outTriangles;
}
