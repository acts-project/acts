// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// STL include(s)
#include <ostream>

// Acts includes
#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

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
class ParametersBase {
 public:

  /// @brief virtual default destructor to allow for inheritance
  virtual ~ParametersBase() = default;

  /// @brief virtual constructor
  virtual ParametersBase* clone() const = 0;

  /// @brief equality operator
  virtual bool operator==(const ParametersBase& rhs) const = 0;

  /// @brief inequality operator
  ///
  /// @return `not (*this == rhs)`
  bool operator!=(const ParametersBase& rhs) const {
    return !(*this == rhs);
  }

  /// @brief access position in global coordinate system
  ///
  /// @return 3D vector with global position
  virtual ActsVectorD<3> position() const = 0;

  /// @brief access momentum in global coordinate system
  ///
  /// @return 3D vector with global momentum
  virtual ActsVectorD<3> momentum() const = 0;

  /// @brief convenience method to retrieve transverse momentum
  double pT() const { return VectorHelpers::perp(momentum()); }

  /// @brief convenience method to retrieve pseudorapidity
  double eta() const { return VectorHelpers::eta(momentum()); }

  /// @brief retrieve electric charge
  ///
  /// @return value of electric charge
  virtual double charge() const = 0;

  /// @brief retrieve time
  ///
  /// @return value of time
  virtual double time() const = 0;

  /// @brief output stream operator
  ///
  /// Prints information about this object to the output stream using the
  /// virtual
  /// TrackParameters::print method.
  ///
  /// @return modified output stream object
  friend std::ostream& operator<<(std::ostream& out,
                                  const ParametersBase& tp) {
    tp.print(out);
    return out;
  }

 protected:
  /// @brief print information to output stream
  ///
  /// @return modified output stream object
  virtual std::ostream& print(std::ostream& sl) const = 0;
};
}  // namespace Acts