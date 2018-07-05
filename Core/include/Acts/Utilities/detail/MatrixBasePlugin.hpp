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
inline Scalar
perp() const
{
  if (this->rows() < 2) return 0.;
  return std::sqrt((*this)[0] * (*this)[0] + (*this)[1] * (*this)[1]);
}

/** perp2 method - perpendicular length squared */
inline Scalar
perp2() const
{
  if (this->rows() < 2) return 0.;
  return ((*this)[0] * (*this)[0] + (*this)[1] * (*this)[1]);
}

inline Scalar
perp2(const MatrixBase<Derived>& vec)
{
  if (this->rows() < 2) return 0.;
  Scalar tot = vec.squaredNorm();
  if (tot > 0) {
    Scalar s = this->dot(vec);
    return this->squaredNorm() - s * s / tot;
  }
  return this->squaredNorm();
}

inline Scalar
perp(const MatrixBase<Derived>& vec)
{
  return std::sqrt(this->perp2(vec));
}

/** phi method */
inline Scalar
phi() const
{
  if (this->rows() < 2) return 0.;
  return std::atan2((*this)[1], (*this)[0]);
}

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
