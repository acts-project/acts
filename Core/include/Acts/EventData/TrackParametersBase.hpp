// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <ostream>

#include "Acts/Utilities/Definitions.hpp"
// Acts includes
#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
// forward declarations
class Surface;

/// @class TrackParametersBase
///
/// @brief base class for track parameters
///
/// This is a base class for neutral and charged track parameters.
/// The position and the momentum are both given in the global coordinate
/// system. The track parameters and their uncertainty are defined in local
/// reference frame which depends on the associated surface
/// of the track parameters.
class TrackParametersBase
{
public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = ActsVector<ParValue_t, Acts::NGlobalPars>;

  /// type of covariance matrix
  using CovMatrix_t = ActsSymMatrix<ParValue_t, Acts::NGlobalPars>;

  /// @brief virtual default destructor to allow for inheritance
  virtual ~TrackParametersBase() = default;

  /// @brief virtual constructor
  virtual TrackParametersBase*
  clone() const = 0;

  /// @brief equality operator
  virtual bool
  operator==(const TrackParametersBase& rhs) const = 0;

  /// @brief inequality operator
  ///
  /// @return `not (*this == rhs)`
  bool
  operator!=(const TrackParametersBase& rhs) const
  {
    return !(*this == rhs);
  }

  /// @brief access position in global coordinate system
  ///
  /// @return 3D vector with global position
  virtual ActsVectorD<3>
  position() const = 0;

  /// @brief access momentum in global coordinate system
  ///
  /// @return 3D vector with global momentum
  virtual ActsVectorD<3>
  momentum() const = 0;

  /// @brief access track parameters
  ///
  /// @return Eigen vector of dimension Acts::NGlobalPars with values of the
  /// track parameters
  ///         (in the order as defined by the ParID_t enumeration)
  ParVector_t
  parameters() const
  {
    return getParameterSet().getParameters();
  }

  /// @brief access track parameter
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter
  ///
  /// @sa ParameterSet::get
  template <ParID_t par>
  ParValue_t
  get() const
  {
    return getParameterSet().template getParameter<par>();
  }

  /// @brief access track parameter uncertainty
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter uncertainty
  template <ParID_t par>
  ParValue_t
  uncertainty() const
  {
    return getParameterSet().template getUncertainty<par>();
  }

  /// @brief access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const CovMatrix_t*
  covariance() const
  {
    return getParameterSet().getCovariance();
  }

  /// @brief convenience method to retrieve transverse momentum
  double
  pT() const
  {
    return VectorHelpers::perp(momentum());
  }

  /// @brief convenience method to retrieve pseudorapidity
  double
  eta() const
  {
    return VectorHelpers::eta(momentum());
  }

  /// @brief retrieve electric charge
  ///
  /// @return value of electric charge
  virtual double
  charge() const = 0;

  /// @brief access associated surface defining the coordinate system for track
  ///        parameters and their covariance
  ///
  /// @return associated surface
  virtual const Surface&
  referenceSurface() const = 0;

  /// @brief access to the measurement frame, i.e. the rotation matrix with
  /// respect to the global coordinate system, in which the local error
  /// is described.
  ///
  /// For planar surface, this is identical to the rotation matrix of the
  /// surface frame, for measurements with respect to a line this has to be
  /// constructed by the point of clostest approach to the line, for
  /// cylindrical surfaces this is (by convention) the tangential plane.
  virtual RotationMatrix3D
  referenceFrame() const = 0;

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual
  /// TrackParameters::print method.
  ///
  /// @return modified output stream object
  friend std::ostream&
  operator<<(std::ostream& out, const TrackParametersBase& tp)
  {
    tp.print(out);
    return out;
  }

  /// @brief access to the internally stored ParameterSet
  ///
  /// @return ParameterSet object holding parameter values and their covariance
  /// matrix
  virtual const FullParameterSet&
  getParameterSet() const = 0;

protected:
  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  virtual std::ostream&
  print(std::ostream& sl) const;
};
}  // namespace Acts