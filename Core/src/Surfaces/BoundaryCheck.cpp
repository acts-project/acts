// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BoundaryCheck.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/BoundaryCheck.hpp"

double Acts::BoundaryCheck::s_cos22 = cos(0.125 * M_PI);
double Acts::BoundaryCheck::s_cos45 = cos(0.25 * M_PI);
double Acts::BoundaryCheck::s_cos67 = cos(0.375 * M_PI);

Acts::BoundaryCheck::BoundaryCheck(bool sCheck)
  : checkLoc0(sCheck)
  , checkLoc1(sCheck)
  , toleranceLoc0(0.)
  , toleranceLoc1(0.)
  , nSigmas(-1)
  , lCovariance(nullptr)
  , bcType(absolute)
{}

Acts::BoundaryCheck::BoundaryCheck(bool   chkL0,
                                   bool   chkL1,
                                   double tloc0,
                                   double tloc1)
  : checkLoc0(chkL0)
  , checkLoc1(chkL1)
  , toleranceLoc0(tloc0)
  , toleranceLoc1(tloc1)
  , nSigmas(-1)
  , lCovariance(nullptr)
  , bcType(absolute)
{
}

Acts::BoundaryCheck::BoundaryCheck(const ActsSymMatrixD<2>& lCov,
                                   double                   nsig,
                                   bool                     chkL0,
                                   bool                     chkL1)
  : checkLoc0(chkL0)
  , checkLoc1(chkL1)
  , toleranceLoc0(0.)
  , toleranceLoc1(0.)
  , nSigmas(nsig)
  , lCovariance(std::make_unique< ActsSymMatrixD<2> >(lCov))
  , bcType(chi2corr)
{
}

Acts::BoundaryCheck::BoundaryCheck(const BoundaryCheck& bCheck)
  : checkLoc0(bCheck.checkLoc0)
  , checkLoc1(bCheck.checkLoc1)
  , toleranceLoc0(bCheck.toleranceLoc0)
  , toleranceLoc1(bCheck.toleranceLoc1)
  , nSigmas(bCheck.nSigmas)
  , lCovariance(nullptr)
  , bcType(bCheck.bcType)
{
  lCovariance = bCheck.lCovariance ? 
      std::make_unique< ActsSymMatrixD<2> >(*bCheck.lCovariance) : nullptr; 
}

Acts::BoundaryCheck& Acts::BoundaryCheck::operator=(const BoundaryCheck& bCheck)
{
  if (this != &bCheck){
      checkLoc0     = bCheck.checkLoc0;
      checkLoc1     = bCheck.checkLoc1;
      toleranceLoc0 = bCheck.toleranceLoc0;
      toleranceLoc1 = bCheck.toleranceLoc1;
      nSigmas       = bCheck.nSigmas;
      lCovariance   = nullptr;
      bcType        = bCheck.bcType;
  }
  return (*this);
}  

