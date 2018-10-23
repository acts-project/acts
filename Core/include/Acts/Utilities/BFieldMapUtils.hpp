// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <vector>
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/Units.hpp"

/// Convenience functions to ease creation of and Acts::InterpolatedBFieldMap
/// and to avoid code duplication. Currently implemented for the two most common
/// formats: rz and xyz.

namespace Acts {

class SolenoidBField;

/// Method to setup the FieldMapper
/// @param localToGlobalBin Function mapping the local bins of r,z to the global
/// bin of the map magnetic field value
///
/// e.g.: we have small grid with the
/// values: r={2,3}, z ={4,5}, the corresponding indices are i (belonging to r)
/// and j (belonging to z), the
/// globalIndex is M (belonging to the value of the magnetic field B(r,z)) and
/// the field map is:
///|   r |    i |    z |    j |   B(r,z) |   M |
///|----:|:----:|:----:|:----:|:--------:|:----|
///|   2 |    0 |    4 |    0 |  2.323   |   0 |
///|   2 |    0 |    5 |    1 |  2.334   |   1 |
///|   3 |    1 |    4 |    0 |  2.325   |   2 |
///|   3 |    1 |    5 |    1 |  2.331   |   3 |
///
/// In this case the function would look like:
/// @code
/// [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
///    return (binsRZ.at(0) * nBinsRZ.at(1) + binsRZ.at(1));
/// }
/// @endcode
/// @param[in] rPos Values of the grid points in r
/// @note The values do not need to be sorted or unique (this will be done
/// inside the function)
/// @param[in] zPos Values of the grid points in z
/// @note The values do not need to be sorted or unique (this will be done
/// inside the function)
/// @param[in] bField The magnetic field values inr r and z for all given grid
/// points stored in a vector
/// @note The function localToGlobalBin determines how the magnetic field was
/// stored in the vector in respect to the grid values
/// @param[in] lengthUnit The unit of the grid points
/// @param[in] BFieldUnit The unit of the magnetic field
/// @param[in] firstQuadrant Flag if set to true indicating that only the first
/// quadrant of the grid points and the BField values has been given and that
/// the BFieldMap should be created symmetrically for all quadrants.
/// e.g. we have the grid values r={0,1} with BFieldValues={2,3} on the r axis.
/// If the flag is set to true the r-axis grid values will be set to {-1,0,1}
/// and the BFieldValues will be set to {3,2,3}.
Acts::InterpolatedBFieldMap::FieldMapper<2, 2>
fieldMapperRZ(const std::function<size_t(std::array<size_t, 2> binsRZ,
                                         std::array<size_t, 2> nBinsRZ)>&
                                          localToGlobalBin,
              std::vector<double>         rPos,
              std::vector<double>         zPos,
              std::vector<Acts::Vector2D> bField,
              double                      lengthUnit    = Acts::units::_mm,
              double                      BFieldUnit    = Acts::units::_T,
              bool                        firstQuadrant = false);

/// Method to setup the FieldMapper
/// @param localToGlobalBin Function mapping the local bins of x,y,z to the
/// global bin of the map magnetic field value
///
/// e.g.: we have small grid with the
/// values: x={2,3}, y={3,4}, z ={4,5}, the corresponding indices are i
/// (belonging to x), j (belonging to y)
/// and k (belonging to z), the globalIndex is M (belonging to the value of the
/// magnetic field B(x,y,z)) and the field map is:
///
///| x   |    i |    y |    j |    z |    k | B(x,y,z) |   M |
///|----:|:----:|:----:|:----:|:----:|:----:|:--------:|:----|
///|   2 |    0 |    3 |    0 |    4 |    0 |  2.323   |   0 |
///|   2 |    0 |    3 |    0 |    5 |    1 |  2.334   |   1 |
///|   2 |    0 |    4 |    1 |    4 |    0 |  2.325   |   2 |
///|   2 |    0 |    4 |    1 |    5 |    1 |  2.331   |   3 |
///|   3 |    1 |    3 |    0 |    4 |    0 |  2.323   |   4 |
///|   3 |    1 |    3 |    0 |    5 |    1 |  2.334   |   5 |
///|   3 |    1 |    4 |    1 |    4 |    0 |  2.325   |   6 |
///|   3 |    1 |    4 |    1 |    5 |    1 |  2.331   |   7 |
///
/// In this case the function would look like:
/// @code
/// [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
///   return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2))
///        + binsXYZ.at(1) * nBinsXYZ.at(2)
///        + binsXYZ.at(2));
/// }
/// @endcode
/// @param[in] xPos Values of the grid points in x
/// @note The values do not need to be sorted or unique (this will be done
/// inside the function)
/// @param[in] yPos Values of the grid points in y
/// @note The values do not need to be sorted or unique (this will be done
/// inside the function)
/// @param[in] zPos Values of the grid points in z
/// @note The values do not need to be sorted or unique (this will be done
/// inside the function)
/// @param[in] bField The magnetic field values inr r and z for all given grid
/// points stored in a vector
/// @note The function localToGlobalBin determines how the magnetic field was
/// stored in the vector in respect to the grid values
/// @param[in] lengthUnit The unit of the grid points
/// @param[in] BFieldUnit The unit of the magnetic field
/// @param[in] firstOctant Flag if set to true indicating that only the first
/// octant of the grid points and the BField values has been given and that
/// the BFieldMap should be created symmetrically for all quadrants.
/// e.g. we have the grid values z={0,1} with BFieldValues={2,3} on the r axis.
/// If the flag is set to true the z-axis grid values will be set to {-1,0,1}
/// and the BFieldValues will be set to {3,2,3}.
Acts::InterpolatedBFieldMap::FieldMapper<3, 3>
fieldMapperXYZ(const std::function<size_t(std::array<size_t, 3> binsXYZ,
                                          std::array<size_t, 3> nBinsXYZ)>&
                                           localToGlobalBin,
               std::vector<double>         xPos,
               std::vector<double>         yPos,
               std::vector<double>         zPos,
               std::vector<Acts::Vector3D> bField,
               double                      lengthUnit  = Acts::units::_mm,
               double                      BFieldUnit  = Acts::units::_T,
               bool                        firstOctant = false);

/// Function which takes an existing SolenoidBField instance and
/// creates a field mapper by sampling grid points from the analytical
/// solenoid field.
///
/// @param rlim pair of r bounds
/// @param zlim pair of z bounds
/// @param nbins pair of bin counts
/// @param field the solenoid field instance
///
/// @return A field mapper instance for use in interpolation.
Acts::InterpolatedBFieldMap::FieldMapper<2, 2>
solenoidFieldMapper(std::pair<double, double> rlim,
                    std::pair<double, double> zlim,
                    std::pair<size_t, size_t> nbins,
                    const SolenoidBField& field);

}  // namespace Acts
