// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <atomic>
#include <vector>

namespace Acts::Experimental {

/// Region-of-interest descriptor with eta/phi/zed bounds.
class RoiDescriptor {
 public:
  /// Iterator type for RoI constituents
  using roi_iterator = std::vector<const RoiDescriptor*>::const_iterator;
  /// convenient
  static constexpr bool FULLSCAN = true;
  /// Flag for region-of-interest mode
  static constexpr bool ROI = false;

  /**
   * @brief constructor
   * @param eta      eta of RoI
   * @param etaMinus eta at rear  of RoI
   * @param etaPlus  eta at front of RoI
   * @param phi      phi of RoI
   * @param phiMinus minimum phi of RoI
   * @param phiPlus  maximum phi of RoI
   * @param zed      zed of RoI
   * @param zedMinus zed at rear  of RoI
   * @param zedPlus  zed at front of RoI
   */
  RoiDescriptor(double eta, double etaMinus, double etaPlus, double phi,
                double phiMinus, double phiPlus, double zed = 0,
                double zedMinus = -s_zedWidthDefault,
                double zedPlus = s_zedWidthDefault);

  // Methods to retrieve data members

  /// Get phi coordinate of RoI center
  /// @return Phi coordinate
  double phi() const { return m_phi; }
  /// Get eta coordinate of RoI center
  /// @return Eta coordinate
  double eta() const { return m_eta; }
  /// Get z coordinate of RoI center
  /// @return Z coordinate
  double zed() const { return m_zed; }

  /// these quantities probably don't need to be used any more
  /// - they are implemented here only because we had them in
  ///   the original legacy interface

  /// Get z at the most forward end of the RoI
  /// @return Z at front
  double zedPlus() const {
    return m_zedPlus;
  }  //!< z at the most forward end of the RoI
  /// Get z at the most backward end of the RoI
  /// @return Z at rear
  double zedMinus() const {
    return m_zedMinus;
  }  //!< z at the most backward end of the RoI

  /// Gets eta at zedPlus
  /// @return Eta at front
  double etaPlus() const { return m_etaPlus; }
  /// Gets eta at zMinus
  /// @return Eta at rear
  double etaMinus() const { return m_etaMinus; }

  /// Gets phiPlus
  /// @return Maximum phi
  double phiPlus() const { return m_phiPlus; }
  /// Gets phiMinus
  /// @return Minimum phi
  double phiMinus() const { return m_phiMinus; }

  /// versioning
  /// Get version identifier
  /// @return Version number
  int version() const { return m_version; }
  /// Set version identifier
  /// @param v Version number
  void version(int v) { m_version = v; }

  /// output
  // virtual operator std::string() const ;

  /// is this a full scan RoI?
  /// @return True if full scan RoI
  bool isFullscan() const { return m_fullscan; }

  /// SuperRoI compatibility methods

  /// am I a SuperRoi?
  /// @return True if composite RoI
  bool composite() const { return m_composite; }
  /// Set composite flag
  /// @param b True to mark as composite
  void setComposite(bool b = true) { m_composite = b; }

  /// Always manage constituents
  /// @return True if managing constituents
  bool manageConstituents() const { return m_manageConstituents; }
  /// Set whether to manage constituents
  /// @param b True to manage constituents
  void manageConstituents(bool b) { m_manageConstituents = b; }

  /// number of constituents
  /// @return Number of constituent RoIs
  unsigned size() const { return m_roiDescriptors.size(); }

  /// find an RoiDescriptor constituent
  /// @param i Index of the constituent
  /// @return The RoI constituent at the given index
  const RoiDescriptor* at(int i) const { return m_roiDescriptors.at(i); }

  /// clear the vector
  void clear() { m_roiDescriptors.clear(); }  // setComposite(false); }

  /// reserve elements in vector
  /// @param s Number of elements to reserve
  void reserve(std::size_t s) { m_roiDescriptors.reserve(s); }

  /// add a RoiDescriptor
  /// @param roi RoI constituent to add
  void push_back(const RoiDescriptor* roi) {
    m_roiDescriptors.push_back(roi);
    setComposite(true);
  }

  /// iterators
  /// @return Begin iterator
  roi_iterator begin() const { return m_roiDescriptors.begin(); }
  /// Get end iterator
  /// @return End iterator
  roi_iterator end() const { return m_roiDescriptors.end(); }

  /// return the gradients
  /// Get dz/dr at the rear of the RoI
  /// @return Gradient dzdr at rear
  double dzdrMinus() const {
    return m_dzdrMinus;
  }  //!<  dz/dr at the rear of the RoI
  /// Get dz/dr at the front of the RoI
  /// @return Gradient dzdr at front
  double dzdrPlus() const { return m_dzdrPlus; }

  /// Get dr/dz at the rear of the RoI
  /// @return Gradient drdz at rear
  double drdzMinus() const { return m_drdzMinus; }
  /// Get dr/dz at the front of the RoI
  /// @return Gradient drdz at front
  double drdzPlus() const { return m_drdzPlus; }

  /// methods to calculate z position at the RoI boundary
  /// at a given radius
  /// Get minimum z at given radius
  /// @param r Radius coordinate
  /// @return Minimum z position
  double zedMin(double r) const;
  /// Get maximum z at given radius
  /// @param r Radius coordinate
  /// @return Maximum z position
  double zedMax(double r) const;

  /// Get z at the most forward end of the RoI at outer radius
  /// @return Z at front outer radius
  double zedOuterPlus() const {
    return m_zedOuterPlus;
  }  //!< z at the most forward end of the RoI
  /// Get z at the most backward end of the RoI at outer radius
  /// @return Z at rear outer radius
  double zedOuterMinus() const {
    return m_zedOuterMinus;
  }  //!< z at the most backward end of the RoI

  /// Get minimum rho at given z
  /// @param z Z coordinate
  /// @return Minimum rho
  double rhoMin(double z) const;
  /// Get maximum rho at given z
  /// @param z Z coordinate
  /// @return Maximum rho
  double rhoMax(double z) const;

  /// Get default z-width
  /// @return Default z width
  static double zedWidthDefault() { return s_zedWidthDefault; }

  /// set default z-width (but only before any RoiDescriptor has been created)
  /// @param d Default z width value
  static void zedWidthDefault(double d);

  // fromn trig
  /// Get RoI identifier
  /// @return RoI ID
  unsigned int roiId() const { return m_roiId; }
  /// Get lvl1 event number
  /// @return Level 1 ID
  unsigned int l1Id() const { return m_l1Id; }
  /// Get lvl1 RoI word
  /// @return RoI word
  unsigned int roiWord() const { return m_roiWord; }

 private:
  /// default parameters - there may be better ways, but this will do
  static std::atomic<double> s_zedWidthDefault;
  /// to ensure default width is only set once at job startup
  static std::atomic<bool> s_firstInstanceCreated;

  float m_phi{};  //!< phi of RoI center
  float m_eta{};  //!< eta of RoI center
  float m_zed{};  //!< zed of RoI center

  float m_phiMinus{};  //!< most negative RoI in azimuthal
  float m_phiPlus{};   //!< most positive RoI in azimuthal
  float m_etaMinus{};  //!< eta of RoI at zedMinus
  float m_etaPlus{};   //!< eta of RoI at zedPlus
  float m_zedMinus{};  //!< z position at most negative position along the
                       //!< beamline
  float
      m_zedPlus{};  //!< z position at most positive position along the beamline

  float m_dzdrMinus{};  //!<  dz/dr at the rear of the RoI
  float m_dzdrPlus{};   //!<  dz/dr at the front of the RoI

  float m_drdzMinus{};  //!<  dr/dz at the rear of the RoI
  float m_drdzPlus{};   //!<  dr/dz at the front of the RoI

  float
      m_zedOuterMinus{};  //!< z at rear of RoI at the outer radius ( = 1100 mm)
  float
      m_zedOuterPlus{};  //!< z at front of RoI at the outer radius ( = 1100 mm)

  bool m_fullscan{};            //!< flag this as a full detector RoI
  bool m_composite{};           //!< flag this as a composite RoI
  bool m_manageConstituents{};  //!< flag to determine whether constituents
                                //!< should be managed

  int m_version{};  //!< transient version identifier

  std::vector<const RoiDescriptor*> m_roiDescriptors;  //!< roi constituents

  // from trig
  unsigned int m_l1Id{0};     //!< lvl1 event number
  unsigned int m_roiId{0};    //!< RoI number
  unsigned int m_roiWord{0};  //!< lvl1 RoI word from which this RoI was
                              //!< initially constructed

  //   std::string str( const RoiDescriptor& d );                           //<!
  //   printing helper std::ostream& operator<<( std::ostream& m, const
  //   RoiDescriptor& d ); //<! printing helper (wraps above)
};

}  // namespace Acts::Experimental
