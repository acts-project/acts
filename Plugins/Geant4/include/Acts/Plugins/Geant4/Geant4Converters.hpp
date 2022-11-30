// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"

#include <memory>
#include <tuple>

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class G4Box;
class G4Trd;
class G4Trap;
class G4Tubs;
class G4VSolid;

namespace Acts {

struct Geant4AlgebraConverter {
  // A potential scalar between Geant4 and ACTS
  ActsScalar scale = 1.;

  /// @brief  Translate a geometry transform: translation only
  ///
  /// @param trans the translation of the Geant4 object
  ///
  /// @return a Acts transform
  Transform3 transform(const G4ThreeVector& trans);

  /// @brief  Translate a geometry transform
  ///
  /// @param rot the rotation of the Geant4 object
  /// @param trans the translation of the Geant4 object
  ///
  /// @return a Acts transform
  Transform3 transform(const G4RotationMatrix& rot, const G4ThreeVector& trans);
};

class AnnulusBounds;
class CylinderBounds;
class RadialBounds;
class RectangleBounds;
class TrapezoidBounds;
class PlanarBounds;

// The following set of converters convert a Geant4 volume shape
// to an ACTS surface bounds object, this is for converting the volume
// based geometry into a surfaced based one.
//
// The obviously smallest expansion/extrusion is condensed to the epsilon
// thin surface.
struct Geant4ShapeConverter {
  /// A scale between Geant4 and ACTS
  ActsScalar scale = 1.;
  /// A description to keep the axis order, if it is set to false
  /// cyclic re-ordering will happen, otherwise axis flip if needed in
  /// order to keep the system right-handed
  bool keepAxisOrder = false;

  /// @brief Convert to cylinder bounds
  ///
  /// @param g4Tubs a Geant4 tube shape
  ///
  /// @return an Acts Cylinder bounds object
  std::shared_ptr<CylinderBounds> cylinderBounds(const G4Tubs& g4Tubs);

  /// @brief Convert to radial bounds
  ///
  /// @param g4Tubs a Geant4 tube shape
  ///
  /// @return an Acts Radial bounds object
  std::shared_ptr<RadialBounds> radialBounds(const G4Tubs& g4Tubs);

  /// @brief Convert to rectangle bounds
  ///
  /// @param g4Box a Geant4 box shape
  ///
  /// @return an ACTS Rectangle bounds shape
  std::tuple<std::shared_ptr<RectangleBounds>, std::array<int, 2u>>
  rectangleBounds(const G4Box& g4Box);

  /// @brief Convert to trapezoid bounds - from Trap
  ///
  /// @param g4Trd a Geant4 trapezoid shape
  ///
  /// @return an ACTS Trapezoid bounds object
  std::tuple<std::shared_ptr<TrapezoidBounds>, std::array<int, 2u>>
  trapezoidBounds(const G4Trd& g4Trd);

  /// @brief Convert to general solid into a planar shape
  ///
  /// @param g4Solid a Geant4 solid shape
  ///
  /// @return an ACTS Planar bounds object
  std::tuple<std::shared_ptr<PlanarBounds>, std::array<int, 2u>> planarBounds(
      const G4VSolid& g4Solid);
};

}  // namespace Acts
