///////////////////////////////////////////////////////////////////
// EnergyLoss.cxx, Acts project
///////////////////////////////////////////////////////////////////

//Fatras module
#include "FatrasUtils/EnergyLoss.h"

// Gaudi
#include "GaudiKernel/MsgStream.h"

#include <iostream>
#include <iomanip>
#include <string>

Acts::EnergyLoss::EnergyLoss() :
  m_deltaE(0.),
  m_sigmaDeltaE(0.),
  m_sigmaMinusDeltaE(0.),
  m_sigmaPlusDeltaE(0.),
  m_mean_ioni(0.),
  m_sig_ioni(0.),
  m_mean_rad(0.),
  m_sig_rad(0.),
  m_length(0.)
{
}

Acts::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                            double sMinusDeltaE, double sPlusDeltaE)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(sMinusDeltaE>0.0?sMinusDeltaE:sigmaDeltaE),
  m_sigmaPlusDeltaE(sPlusDeltaE>0.0?sPlusDeltaE:sigmaDeltaE),
  m_mean_ioni(0.),
  m_sig_ioni(0.),
  m_mean_rad(0.),
  m_sig_rad(0.),
  m_length(-std::numeric_limits<double>::min())
{ }

Acts::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                            double sMinusDeltaE, double sPlusDeltaE, 
                            double mean_ioni, double sig_ioni,
                double mean_rad, double sig_rad, double length)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(sMinusDeltaE),
  m_sigmaPlusDeltaE(sPlusDeltaE),
  m_mean_ioni(mean_ioni),
  m_sig_ioni(sig_ioni),
  m_mean_rad(mean_rad),
  m_sig_rad(sig_rad),
  m_length(length)
{ }

Acts::EnergyLoss::EnergyLoss(double deltaE, double sigmaDeltaE,
                            double mean_ioni, double sig_ioni,
                double mean_rad, double sig_rad)  :
  m_deltaE(deltaE),
  m_sigmaDeltaE(sigmaDeltaE),
  m_sigmaMinusDeltaE(0.),
  m_sigmaPlusDeltaE(0.),
  m_mean_ioni(mean_ioni),
  m_sig_ioni(sig_ioni),
  m_mean_rad(mean_rad),
  m_sig_rad(sig_rad),
  m_length(-std::numeric_limits<double>::min())
{ }

// Overload of << operator for MsgStream and std::ostream for debug output 
MsgStream& Acts::operator << ( MsgStream& sl, const Acts::EnergyLoss& eloss)
{ 
  sl << "EnergyLoss :   ( delta(E), sigma(dE) ) = \t"
     << "("<< eloss.deltaE()<<", \t"<< eloss.sigmaDeltaE() << ")";
  return sl;
}

std::ostream& Acts::operator << ( std::ostream& sl, const Acts::EnergyLoss& eloss)
{ 
  sl << "EnergyLoss :   ( delta(E), sigma(dE) ) = \t"
     << "("<< eloss.deltaE()<<", \t"<< eloss.sigmaDeltaE() << ")";
  return sl; 
}
