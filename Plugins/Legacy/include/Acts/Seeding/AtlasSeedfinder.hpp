// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AtlasSeedfinder.hpp Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Seeding/LegacyInternalSeed.hpp"
#include "Acts/Seeding/SPForSeed.hpp"

#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
namespace Legacy {
template <typename SpacePoint>
class AtlasSeedfinder {
  struct Config {
    // UNIT AS RETURNED BY m_fieldService->getField() default value in ATLAS
    // was
    // 5. Unit is kilo-Tesla
    //  double bFieldInZ = 5.;
    double bFieldInZ = 0.00208;

    double SCT_rMin = 200.;

    double beamPosX = 0;
    double beamPosY = 0;
    double beamPosZ = 0;
    double beamTiltX = 0;
    double beamTiltY = 0;
  };
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

 public:
  ///////////////////////////////////////////////////////////////////
  // Standard tool methods
  ///////////////////////////////////////////////////////////////////

  AtlasSeedfinder();
  virtual ~AtlasSeedfinder();

  ///////////////////////////////////////////////////////////////////
  // Methods to initialize tool for new event or region
  ///////////////////////////////////////////////////////////////////

  template <class RandIter>
  void newEvent(int, RandIter, RandIter);

  //////////////////////////////////////////////////////////////////
  // Method to initialize seeds production
  //////////////////////////////////////////////////////////////////
  void find3Sp();

  ///////////////////////////////////////////////////////////////////
  // Iterator through seeds pseudo collection produced accordingly
  // methods find
  ///////////////////////////////////////////////////////////////////

  const Seed<SpacePoint>* next();

  ///////////////////////////////////////////////////////////////////
  // Configuration
  ///////////////////////////////////////////////////////////////////

  const Config m_config;

 protected:
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  AtlasSeedfinder(const AtlasSeedfinder<SpacePoint>&) = delete;
  AtlasSeedfinder<SpacePoint>& operator=(const AtlasSeedfinder<SpacePoint>&) =
      delete;
  //@}
  ///////////////////////////////////////////////////////////////////
  // Protected data and methods
  ///////////////////////////////////////////////////////////////////

  bool m_endlist;
  bool m_checketa;
  bool m_isvertex;
  int m_nprint;
  int m_nlist;
  int m_maxsize;
  int m_state;
  // event number since tool init
  int m_iteration;
  float m_etamin, m_etamax;
  float m_drmin, m_drminv;
  float m_drmax;
  float m_dzdrmin0;
  float m_dzdrmax0;
  float m_dzdrmin;
  float m_dzdrmax;
  float m_zmin;
  float m_zmax;
  float m_zminU;
  float m_zmaxU;
  float m_zminB;
  float m_zmaxB;
  float m_ftrig;
  float m_ftrigW;
  // maximum radius of outermost detector element
  float r_rmax;
  // size of one r-slice
  float r_rstep;

  float m_dzver;
  float m_dzdrver;
  float m_diver;
  float m_diverpps;
  float m_diversss;
  float m_divermax;
  float m_dazmax;
  float m_ptmin;
  float m_ipt;
  float m_ipt2;
  float m_COF;
  float m_K;
  float m_ipt2K;
  float m_ipt2C;
  float m_COFK;
  float m_umax;
  // number of r-slices
  int r_size;
  int r_first;
  int rf_size;
  int rfz_size;
  std::list<SPForSeed<SpacePoint>*>* r_Sorted;
  std::list<SPForSeed<SpacePoint>*> rfz_Sorted[583];
  std::list<SPForSeed<SpacePoint>*> l_spforseed;
  typename std::list<SPForSeed<SpacePoint>*>::iterator i_spforseed;
  typename std::list<SPForSeed<SpacePoint>*>::iterator m_rMin;

  int m_nsaz, m_nsazv;
  int m_fNmax, m_fvNmax;
  int m_fNmin, m_fvNmin;
  int m_zMin;
  // m_nr: number of bins used in r_Sorted; r_index: index of all used bins in
  // r_Sorted; r_map is number of SP in each bin in r_Sorted
  int m_nr;
  int* r_index;
  int* r_map;
  int m_nrfz, rfz_index[583], rfz_map[583];
  int rfz_b[583], rfz_t[593], rfz_ib[583][9], rfz_it[583][9];
  float m_sF;

  ///////////////////////////////////////////////////////////////////
  // Tables for 3 space points seeds search
  ///////////////////////////////////////////////////////////////////

  int m_maxsizeSP;
  SPForSeed<SpacePoint>** m_SP;
  float* m_Zo;
  float* m_Tz;
  float* m_R;
  float* m_U;
  float* m_V;
  float* m_Er;

  Seed<SpacePoint>* m_seedOutput;

  std::list<InternalSeed<SpacePoint>*> l_seeds;
  typename std::list<InternalSeed<SpacePoint>*>::iterator i_seed;
  typename std::list<InternalSeed<SpacePoint>*>::iterator i_seede;

  std::multimap<float, InternalSeed<SpacePoint>*> m_seeds;
  typename std::multimap<float, InternalSeed<SpacePoint>*>::iterator m_seed;

  std::multimap<float, InternalSeed<SpacePoint>*> m_mapOneSeeds;
  InternalSeed<SpacePoint>* m_OneSeeds;
  int m_maxOneSize;
  int m_nOneSeeds;
  int m_fillOneSeeds;
  std::vector<std::pair<float, SPForSeed<SpacePoint>*>> m_CmSp;

  ///////////////////////////////////////////////////////////////////
  // Beam geometry
  ///////////////////////////////////////////////////////////////////

  float m_xbeam;  // x-center of beam-axis
  float m_ybeam;  // y-center of beam-axis
  float m_zbeam;  // z-center of beam-axis

  ///////////////////////////////////////////////////////////////////
  // Protected methods
  ///////////////////////////////////////////////////////////////////

  void buildFrameWork();
  void buildBeamFrameWork();

  SPForSeed<SpacePoint>* newSpacePoint(SpacePoint* const&);

  void newOneSeed(SPForSeed<SpacePoint>*&, SPForSeed<SpacePoint>*&,
                  SPForSeed<SpacePoint>*&, float, float);

  void newOneSeedWithCurvaturesComparison(SPForSeed<SpacePoint>*&,
                                          SPForSeed<SpacePoint>*&, float);

  void fillSeeds();
  void fillLists();
  void erase();
  void production3Sp();
  void production3Sp(typename std::list<SPForSeed<SpacePoint>*>::iterator*,
                     typename std::list<SPForSeed<SpacePoint>*>::iterator*,
                     typename std::list<SPForSeed<SpacePoint>*>::iterator*,
                     typename std::list<SPForSeed<SpacePoint>*>::iterator*, int,
                     int, int&);

  void findNext();
  bool isZCompatible(float&);
  void convertToBeamFrameWork(SpacePoint* const&, float*);
};

///////////////////////////////////////////////////////////////////
// Inline methods
///////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline const Seed<SpacePoint>* AtlasSeedfinder<SpacePoint>::next() {
  do {
    if (i_seed == i_seede) {
      findNext();
      if (i_seed == i_seede) {
        return 0;
      }
    }
    ++i_seed;
  } while (!(*m_seed++).second->set3(*m_seedOutput));
  return (m_seedOutput);
}

template <typename SpacePoint>
inline bool AtlasSeedfinder<SpacePoint>::isZCompatible(float& Zv) {
  if (Zv < m_zminU || Zv > m_zmaxU) {
    return false;
  } else {
    return true;
  }
}

///////////////////////////////////////////////////////////////////
// New space point for seeds
///////////////////////////////////////////////////////////////////

template <typename SpacePoint>
inline SPForSeed<SpacePoint>* AtlasSeedfinder<SpacePoint>::newSpacePoint(
    SpacePoint* const& sp) {
  SPForSeed<SpacePoint>* sps;

  float r[3];
  convertToBeamFrameWork(sp, r);

  if (m_checketa) {
    // filter SP outside of eta-range
    float z = (fabs(r[2]) + m_zmax);
    float x = r[0] * m_dzdrmin;
    float y = r[1] * m_dzdrmin;
    if ((z * z) < (x * x + y * y)) {
      return 0;
    }
  }

  if (i_spforseed != l_spforseed.end()) {
    sps = (*i_spforseed++);
    sps->set(sp, r);
  } else {
    l_spforseed.push_back((sps = new SPForSeed<SpacePoint>(sp, r)));
    i_spforseed = l_spforseed.end();
  }

  return sps;
}

///////////////////////////////////////////////////////////////////
// Object-function for curvature seeds comparison
///////////////////////////////////////////////////////////////////

class comCurvature {
 public:
  template <typename SpacePoint>
  bool operator()(
      const std::pair<float, Acts::Legacy::SPForSeed<SpacePoint>*>& i1,
      const std::pair<float, Acts::Legacy::SPForSeed<SpacePoint>*>& i2) {
    return i1.first < i2.first;
  }
};
}  // namespace Legacy
}  // namespace Acts
#include "Acts/Seeding/AtlasSeedfinder.ipp"
