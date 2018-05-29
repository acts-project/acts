// This file is part of the Acts project.
//
// Copyright (C) 2016 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TransportJacobian.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/EventData/TransportJacobian.hpp"
#include <iomanip>
#include <iostream>

///////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////

Acts::TransportJacobian::TransportJacobian(const double* J)
  : ActsMatrixD<5, 5>()
{
  (*this)(0, 0) = J[0];
  (*this)(0, 1) = J[1];
  (*this)(0, 2) = J[2];
  (*this)(0, 3) = J[3];
  (*this)(0, 4) = J[4];
  (*this)(1, 0) = J[5];
  (*this)(1, 1) = J[6];
  (*this)(1, 2) = J[7];
  (*this)(1, 3) = J[8];
  (*this)(1, 4) = J[9];
  (*this)(2, 0) = J[10];
  (*this)(2, 1) = J[11];
  (*this)(2, 2) = J[12];
  (*this)(2, 3) = J[13];
  (*this)(2, 4) = J[14];
  (*this)(3, 0) = J[15];
  (*this)(3, 1) = J[16];
  (*this)(3, 2) = J[17];
  (*this)(3, 3) = J[18];
  (*this)(3, 4) = J[19];
  (*this)(4, 0) = J[20];
  (*this)(4, 1) = J[21];
  (*this)(4, 2) = J[22];
  (*this)(4, 3) = J[23];
  (*this)(4, 4) = J[24];
}

Acts::TransportJacobian::TransportJacobian(const ActsMatrixD<5, 5>& J)
  : ActsMatrixD<5, 5>(J)
{
}

///////////////////////////////////////////////////////////////////
// Dumps relevant information into the ostream
///////////////////////////////////////////////////////////////////

std::ostream&
Acts::operator<<(std::ostream& sl, const Acts::TransportJacobian& J)
{
  sl << "|-------------|-------------|-------------|-------------|-------------"
        "|-------------|"
     << std::endl;
  sl << "| Jacobian(A) | Old   /dL1  |       /dL2  |      /dPhi  |      /dThe  "
        "|       /dCM  |"
     << std::endl;
  sl << "|-------------|-------------|-------------|-------------|-------------"
        "|-------------|"
     << std::endl;

  sl << "|  New  dL1 / |" << std::setw(12) << std::setprecision(5) << J(0, 0)
     << " |" << std::setw(12) << std::setprecision(5) << J(0, 1) << " |"
     << std::setw(12) << std::setprecision(5) << J(0, 2) << " |"
     << std::setw(12) << std::setprecision(5) << J(0, 3) << " |"
     << std::setw(12) << std::setprecision(5) << J(0, 4) << " |" << std::endl;
  sl << "|       dL2 / |" << std::setw(12) << std::setprecision(5) << J(1, 0)
     << " |" << std::setw(12) << std::setprecision(5) << J(1, 1) << " |"
     << std::setw(12) << std::setprecision(5) << J(1, 2) << " |"
     << std::setw(12) << std::setprecision(5) << J(1, 3) << " |"
     << std::setw(12) << std::setprecision(5) << J(1, 4) << " |" << std::endl;
  sl << "|       dPhi/ |" << std::setw(12) << std::setprecision(5) << J(2, 0)
     << " |" << std::setw(12) << std::setprecision(5) << J(2, 1) << " |"
     << std::setw(12) << std::setprecision(5) << J(2, 2) << " |"
     << std::setw(12) << std::setprecision(5) << J(2, 3) << " |"
     << std::setw(12) << std::setprecision(5) << J(2, 4) << " |" << std::endl;
  sl << "|       dThe/ |" << std::setw(12) << std::setprecision(5) << J(3, 0)
     << " |" << std::setw(12) << std::setprecision(5) << J(3, 1) << " |"
     << std::setw(12) << std::setprecision(5) << J(3, 2) << " |"
     << std::setw(12) << std::setprecision(5) << J(3, 3) << " |"
     << std::setw(12) << std::setprecision(5) << J(3, 4) << " |" << std::endl;
  sl << "|       dCM / |"

     << std::setw(12) << std::setprecision(5) << J(4, 0) << " |"
     << std::setw(12) << std::setprecision(5) << J(4, 1) << " |"
     << std::setw(12) << std::setprecision(5) << J(4, 2) << " |"
     << std::setw(12) << std::setprecision(5) << J(4, 3) << " |"
     << std::setw(12) << std::setprecision(5) << J(4, 4) << " |" << std::endl;
  sl << "|-------------|-------------|-------------|-------------|-------------"
        "|-------------|"
     << std::endl;
  return sl;
}
