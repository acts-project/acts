// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// BinUtility.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <array>
#include <iostream>
#include <memory>
#include <vector>
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

///  @class BinUtility
///
/// The BinUtility class that translated global and local position into a bins
/// of a BinnedArray, most performant is equidistant binning without a
/// transform,
/// however, optionally a transform can be provided, e.g. for binning on shifted
/// object, the transform is usually shared with the geometric object the Array
/// is
/// defined on, for performance reasons, also the inverse transform is stored.
///
class BinUtility {
 public:
  /// Constructor for equidistant
  BinUtility() : m_binningData(), m_transform(nullptr), m_itransform(nullptr) {
    m_binningData.reserve(3);
  }

  /// Constructor from BinningData directly
  ///
  /// @param bData is the provided binning data
  /// @param tForm is the (optional) transform
  BinUtility(const BinningData& bData,
             const std::shared_ptr<const Transform3D>& tForm = nullptr)
      : m_binningData(),
        m_transform(tForm),
        m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr) {
    m_binningData.reserve(3);
    m_binningData.push_back(bData);
  }

  /// Constructor for equidistant
  ///
  /// @param bins is the number of bins
  /// @param min in the minimal value
  /// @param max is the maximal value
  /// @param opt is the binning option : open, closed
  /// @param value is the binninb value : binX, binY, binZ, etc.
  /// @param tForm is the (optional) transform
  BinUtility(size_t bins, float min, float max, BinningOption opt = open,
             BinningValue value = binX,
             const std::shared_ptr<const Transform3D>& tForm = nullptr)
      : m_binningData(),
        m_transform(tForm),
        m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr) {
    m_binningData.reserve(3);
    m_binningData.push_back(BinningData(opt, value, bins, min, max));
  }

  /// Constructor for arbitrary
  ///
  /// @param bValues is the boundary values of the binning
  /// @param opt is the binning option : open, closed
  /// @param value is the binninb value : binX, binY, binZ, etc.
  /// @param tForm is the (optional) transform
  BinUtility(std::vector<float>& bValues, BinningOption opt = open,
             BinningValue value = binPhi,
             const std::shared_ptr<const Transform3D>& tForm = nullptr)
      : m_binningData(),
        m_transform(tForm),
        m_itransform(tForm ? new Transform3D(tForm->inverse()) : nullptr) {
    m_binningData.reserve(3);
    m_binningData.push_back(BinningData(opt, value, bValues));
  }

  /// Copy constructor
  ///
  /// @param sbu is the source bin utility
  BinUtility(const BinUtility& sbu)
      : m_binningData(sbu.m_binningData),
        m_transform(sbu.m_transform),
        m_itransform(sbu.m_transform
                         ? new Transform3D(sbu.m_transform->inverse())
                         : nullptr) {}

  /// Assignment operator
  ///
  /// @param sbu is the source bin utility
  BinUtility& operator=(const BinUtility& sbu) {
    if (this != &sbu) {
      m_binningData = sbu.m_binningData;
      m_transform = sbu.m_transform;
      m_itransform = sbu.m_transform
                         ? std::unique_ptr<const Transform3D>(
                               new Transform3D(sbu.m_transform->inverse()))
                         : nullptr;
    }
    return (*this);
  }

  /// Operator++ to make multidimensional BinUtility
  ///
  /// @param gbu is the additional BinUtility to be chosen
  BinUtility& operator+=(const BinUtility& gbu) {
    const std::vector<BinningData>& bData = gbu.binningData();

    if (m_transform == nullptr && gbu.transform() != nullptr) {
      // use other transform
      m_transform = gbu.transform();
      m_itransform =
          std::make_unique<const Transform3D>(m_transform->inverse());
    } else if (m_transform != nullptr && gbu.transform() != nullptr) {
      // combine two existing transform
      // note that this might lead to undesired behaviour of the combined
      // BinUtility
      m_transform = std::make_shared<const Transform3D>((*m_transform) *
                                                        (*gbu.transform()));
      m_itransform =
          std::make_unique<const Transform3D>(m_transform->inverse());
    }  // else {
    // only this BU has transform, just keep it.
    //}

    if (m_binningData.size() + bData.size() > 3) {
      throw "BinUtility does not support dim > 3";
    }
    m_binningData.insert(m_binningData.end(), bData.begin(), bData.end());
    return (*this);
  }

  /// Virtual Destructor
  ~BinUtility() = default;

  /// return the binning data vector
  const std::vector<BinningData>& binningData() const { return m_binningData; }

  /// Return the total number of bins
  size_t bins() const { return bins(0) * bins(1) * bins(2); }

  /// Bin-triple fast access
  ///
  /// - calculate the bin triple with one transform
  ///
  /// @param position is the 3D position to be evaluated
  ///
  /// @return is the bin value in 3D
  std::array<size_t, 3> binTriple(const Vector3D& position) const {
    /// transform or not
    const Vector3D& bPosition =
        m_itransform ? Vector3D((*m_itransform) * position) : position;
    // get the dimension
    size_t mdim = m_binningData.size();
    /// now get the bins
    size_t bin0 = m_binningData[0].searchGlobal(bPosition);
    size_t bin1 = mdim > 1 ? m_binningData[1].searchGlobal(bPosition) : 0;
    size_t bin2 = mdim > 2 ? m_binningData[2].searchGlobal(bPosition) : 0;
    /// return the triple
    return {{bin0, bin1, bin2}};
  }

  /// Bin from a 3D vector (already in binning frame)
  /// @todo - optionally the itransform is applied
  ///
  /// @param position is the 3D position to be evaluated
  /// @param ba is the bin dimension
  ///
  /// @return is the bin value
  size_t bin(const Vector3D& position, size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      return 0;
    }
    size_t bEval =
        m_itransform
            ? m_binningData[ba].searchGlobal((*m_itransform) * position)
            : m_binningData[ba].searchGlobal(position);
    return bEval;
  }

  /// Bin neighbour range
  /// this method calls the increment/decreement methods
  /// the bin itself is also contained, so if not an edge-case
  /// this would be
  ///  | n | c | p |
  ///
  /// @param position is the position for the neighbour Range test
  /// @param ba is the binning accessor
  ///
  /// @return a vector of neighbour sizes
  std::vector<size_t> neighbourRange(const Vector3D& position,
                                     size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      return {0};
    }
    std::vector<size_t> neighbourRange;
    size_t cbin = bin(position, ba);
    size_t pbin = cbin;
    size_t nbin = cbin;
    if (m_binningData[ba].decrement(pbin)) {
      neighbourRange.push_back(pbin);
    }
    neighbourRange.push_back(cbin);
    if (m_binningData[ba].increment(nbin) && nbin != pbin) {
      neighbourRange.push_back(nbin);
    }
    return neighbourRange;
  }

  /// Return the oder direction for fast interlinking
  ///
  /// @param position is the global position for the next search
  /// @param direction is the global position for the next search
  /// @param ba is the bin accessor
  ///
  /// @todo the
  ///
  /// @return the next bin
  int nextDirection(const Vector3D& position, const Vector3D& direction,
                    size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      return 0;
    }
    return m_binningData[ba].nextDirection(position, direction);
  }

  /// Bin from a 2D vector (following local parameters defintitions)
  /// - no optional transform applied
  /// - USE WITH CARE !!
  ///
  /// You need to make sure that the local position is actually in the binning
  /// frame of the BinUtility
  ///
  /// @param lposition is the local position to be set
  /// @param ba is the bin dimension
  ///
  /// @return bin calculated from local
  size_t bin(const Vector2D& lposition, size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      return 0;
    }
    return m_binningData[ba].searchLocal(lposition);
  }
  /// Check if bin is inside from Vector2D - optional transform applied
  ///
  /// @param position is the global position to be evaluated
  /// @return is a boolean check
  bool inside(const Vector3D& position) const {
    /// transform or not
    const Vector3D& bPosition =
        m_itransform ? Vector3D((*m_itransform) * position) : position;
    // loop and break
    for (auto& bData : m_binningData) {
      if (!(bData.inside(bPosition))) {
        return false;
      }
    }
    // survived all the checks
    return true;
  }

  /// Check if bin is inside from Vector2D - no optional transform applied
  ///
  /// @param lposition is the local position to be evaluated
  /// @return is a boolean check
  bool inside(const Vector2D& lposition) const {
    return true;
    std::vector<BinningData>::const_iterator bdIter = m_binningData.begin();
    for (; bdIter != m_binningData.end(); ++bdIter) {
      if (!(*bdIter).inside(lposition)) {
        return false;
      }
    }
    return true;
  }

  /// First bin maximal value
  /// @return the dimenstion of the binning data
  size_t dimensions() const { return m_binningData.size(); }

  /// First bin maximal value
  ///
  /// @param ba is the binaccessor
  ///
  /// @return size_t is the maximal bin of the accessor entry
  size_t max(size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      return 0;
    }
    return (m_binningData[ba].bins() - 1);
  }

  /// Number of bins
  ///
  /// @param ba is the binaccessor
  ///
  /// @return size_t is the bins of the accessor entry
  size_t bins(size_t ba) const {
    if (ba >= m_binningData.size()) {
      return 1;
    }
    return (m_binningData[ba].bins());
  }

  /// Transform applied to global positions before lookup
  ///
  /// @return Shared pointer to transform
  std::shared_ptr<const Transform3D> transform() const { return m_transform; }

  /// The type/value of the binning
  ///
  /// @param ba is the binaccessor
  ///
  /// @return the binning value of the accessor entry
  BinningValue binningValue(size_t ba = 0) const {
    if (ba >= m_binningData.size()) {
      throw "dimension out of bounds";
    }
    return (m_binningData[ba].binvalue);
  }

  /// Serialize the bin triple
  /// - this creates a simple size_t from a triple object
  ///
  /// @param bin is the bin to be serialized
  size_t serialize(const std::array<size_t, 3>& bin) const {
    size_t serializedBin = bin[0];
    if (m_binningData.size() == 2) {
      serializedBin += bin[1] * m_binningData[0].bins();
    } else if (m_binningData.size() == 3) {
      serializedBin +=
          (bin[1] * m_binningData[0].bins() * bin[2] * m_binningData[1].bins());
    }
    return serializedBin;
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const {
    sl << "BinUtility for " << m_binningData.size()
       << "- dimensional array:" << std::endl;
    std::vector<BinningData>::const_iterator bdIter = m_binningData.begin();
    for (size_t ibd = 0; bdIter != m_binningData.end(); ++bdIter, ++ibd) {
      sl << "dimension     : " << ibd << std::endl;
      sl << " - type       : " << size_t((*bdIter).type) << std::endl;
      sl << " - option     : " << size_t((*bdIter).option) << std::endl;
      sl << " - value      : " << size_t((*bdIter).binvalue) << std::endl;
      sl << " - bins       : " << (*bdIter).bins() << std::endl;
      sl << " - min/max    : " << (*bdIter).min << " / " << (*bdIter).max
         << std::endl;
      if ((*bdIter).type == equidistant) {
        sl << " - step       : " << (*bdIter).step << std::endl;
      }
      sl << " - boundaries : | ";
      std::vector<float>::const_iterator bIter = (*bdIter).boundaries().begin();
      for (; bIter != (*bdIter).boundaries().end(); ++bIter) {
        sl << (*bIter) << " | ";
      }
      sl << std::endl;
    }
    return sl;
  }

 private:
  std::vector<BinningData> m_binningData;           /// vector of BinningData
  std::shared_ptr<const Transform3D> m_transform;   /// shared transform
  std::unique_ptr<const Transform3D> m_itransform;  /// unique inverse transform
};

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const BinUtility& bgen);

}  // namespace Acts
