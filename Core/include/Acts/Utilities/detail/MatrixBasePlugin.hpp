// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AmgMatrixPlugin.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

// ------- Methods for 3D vector type objects ---------------------- //

/** perp method - perpendicular length */

/** perp2 method - perpendicular length squared */



/** phi method */

/** theta method */
inline Scalar
theta() const
{
  if (this->rows() < 3) return 0.;
  return std::atan2(
      std::sqrt((*this)[0] * (*this)[0] + (*this)[1] * (*this)[1]), (*this)[2]);
}

/** pseudorapidity  method */
inline Scalar
eta() const
{
  return -std::log(std::tan(this->theta() * .5));  // TODO: slow
}
