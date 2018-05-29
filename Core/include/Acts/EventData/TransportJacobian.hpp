// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TransportJacobian.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class TransportJacobian
///
///   This class represents the jacobian for transforming initial track
/// parameters
///   to new track parameters after propagation to new surface.
///   Initial track parameters:           loc0  ,loc1  ,Phi  ,Theta  ,qp
///   Track parameters after propagation: loc0_n,loc1_n,Phi_n,Theta_n,qp_n
///
///   Jacobian is matrix (5x5) with derivatives
///
///          0                1                2               3 4
///   0 d loc0_n /d loc0 d loc0_n /d loc1 d loc0_n /d Phi d loc0_n /d Theta d
/// loc0_n /d qp
///   1 d loc1_n /d loc0 d loc1_n /d loc1 d loc1_n /d Phi d loc1_n /d Theta d
/// loc1_n /d qp
///   2 d Phi_n  /d loc0 d Phi_n  /d loc1 d Phi_n  /d Phi d Phi_n  /d Theta d
/// Phi_n  /d qp
///   3 d Theta_n/d loc0 d Theta_n/d loc1 d Theta_n/d Phi d Theta_n/d Theta d
/// Theta_n/d qp
///   4 d qp_n   /d loc0 d qp_n   /d loc1 d qp_n   /d Phi d qp   _n/d Theta d
/// qp_n   /d qp
///
///   ^ ---> second index
///   |
///   | first index
///
///
class TransportJacobian : public ActsMatrixD<5, 5>
{
public:
  /// Constructor
  ///
  /// @param tdata is the double array for containing the jacobian entries
  TransportJacobian(const double* tdata);
  /// Constructor
  ///
  /// @param tdata is the matrix containing the jacobian entries
  TransportJacobian(const ActsMatrixD<5, 5>& tdata);

  /// Destructor
  virtual ~TransportJacobian(){};
};

/// Overload of << operator for std::ostream
std::ostream&
operator<<(std::ostream& sl, const TransportJacobian& jac);

}  // end of namespace