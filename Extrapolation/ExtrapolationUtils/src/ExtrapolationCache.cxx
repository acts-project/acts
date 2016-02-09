///////////////////////////////////////////////////////////////////
// ExtrapolationCache.cxx 
///////////////////////////////////////////////////////////////////

#include "ExtrapolationUtils/ExtrapolationCache.h" 
#include "MaterialOnTrack/EnergyLoss.h" 


Ats::ExtrapolationCache::ExtrapolationCache() : 
   m_x0tot(0), m_eloss(0) { }

Ats::ExtrapolationCache::ExtrapolationCache(double x0tot) : 
   m_x0tot(x0tot), m_eloss(0) { }

Ats::ExtrapolationCache::ExtrapolationCache(double x0tot, EnergyLoss* eloss): 
   m_x0tot(x0tot),
   m_eloss(eloss) { }

Ats::ExtrapolationCache::ExtrapolationCache(const Ats::ExtrapolationCache& cache) : 
   m_x0tot(cache.m_x0tot), 
   m_eloss(cache.m_eloss) { }


