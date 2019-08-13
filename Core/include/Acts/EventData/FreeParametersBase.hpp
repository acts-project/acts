// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <ostream>

// Acts includes
#include "Acts/EventData/ParametersBase.hpp"
#include "Acts/EventData/FreeParameterSet.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class TrackParametersBase
///
/// @brief base class for track parameters
///
/// This is a base class for neutral and charged track parameters.
/// The position and the momentum are both given in the global coordinate
/// system. The track parameters and their uncertainty are defined in local
/// reference frame which depends on the associated surface
/// of the track parameters.
class FreeParametersBase : public ParametersBase {
 public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = FreeVector;

  /// type of covariance matrix
  using CovMatrix_t = FreeSymMatrix;

  /// @brief access track parameters
  ///
  /// @return Eigen vector of dimension Acts::BoundParsDim with values of the
  /// track parameters
  ///         (in the order as defined by the ParID_t enumeration)
  ParVector_t parameters() const = 0;

  /// @brief access track parameter
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter
  ///
  /// @sa ParameterSet::get
  template <FreeID_t par>
  ParValue_t get() const  = 0;

  /// @brief access track parameter uncertainty
  ///
  /// @tparam par identifier of track parameter which is to be retrieved
  ///
  /// @return value of the requested track parameter uncertainty
  template <FreeID_t par>
  ParValue_t uncertainty() const = 0;
  
  /// @brief access covariance matrix of track parameters
  ///
  /// @note The ownership of the covariance matrix is @b not transferred with
  /// this call.
  ///
  /// @return raw pointer to covariance matrix (can be a nullptr)
  ///
  /// @sa ParameterSet::getCovariance
  const CovMatrix_t* covariance() const = 0;
  
protected:
  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  virtual std::ostream& print(std::ostream& sl) const;
};
}  // namespace Acts