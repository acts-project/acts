// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/// @note This file is foreseen for the `Geometry` module to replace `Extent`

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Enumerate.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Range1D.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <iostream>
#include <vector>

namespace Acts {

using Envelope = std::array<ActsScalar, 2>;
using ExtentEnvelope = std::array<std::array<ActsScalar, 2>, binValues>;

static Envelope zeroEnvelope = {0., 0};
static ExtentEnvelope zeroEnvelopes = {zeroEnvelope, zeroEnvelope, zeroEnvelope,
                                       zeroEnvelope, zeroEnvelope, zeroEnvelope,
                                       zeroEnvelope, zeroEnvelope};

class GeometricExtent {
 public:
  /// Constructor with (optional) @param envelope
  GeometricExtent(const ExtentEnvelope& envelope = zeroEnvelopes);

  /// Copy constructor
  GeometricExtent(const GeometricExtent& rhs) = default;

  /// Extend with a vector
  ///
  /// @param vtx the vertex to be used for extending
  /// @param bValues the binning values
  /// @param applyEnv boolean to steer if envelope should be applied
  void extend(const Vector3& vtx,
              const std::vector<BinningValue>& bValues = s_binningValues,
              bool applyEnv = true);

  /// Extend with a set of vectors by iterators
  ///
  /// @param start the start iterator of the loop
  /// @param end the end iterator of the loop
  /// @param bValues the binning values
  /// @param applyEnv boolean to steer if envelope should be applied
  template <typename vector_iterator_t>
  void extend(const vector_iterator_t& start, const vector_iterator_t& end,
              const std::vector<BinningValue>& bValues = s_binningValues,
              bool applyEnv = true) {
    for (vector_iterator_t vIt = start; vIt < end; ++vIt) {
      extend(*vIt, bValues, applyEnv);
    }
  }

  /// Extend with another geometric extent, usually pushes the
  /// current range to the boundaries of the @param rhs extent,
  /// unless the current extent is already bigger.
  ///
  /// @note the @param rhs extent can also simply set an envelope
  /// which then is applied to the current one
  ///
  /// @param rhs is the other source Extent
  /// @param bValues the binning values
  /// @param applyEnv boolean to steer if envelope should be applied
  ///        on the constraint values, if only an envelope is given
  ///        but the value not constraint, then it is always applied
  void extend(const GeometricExtent& rhs,
              const std::vector<BinningValue>& bValues = s_binningValues,
              bool applyEnv = true);

  /// Set a range
  ///
  /// @param bValue the binning identification
  /// @param min the minimum parameter
  /// @param max the maximum parameter
  void set(BinningValue bValue, ActsScalar min, ActsScalar max);

  /// (re-)Set the envelope
  ///
  /// @param the new envelope to be set
  void setEnvelope(const ExtentEnvelope& envelope = zeroEnvelopes);

  /// Return the individual 1-dimensional range
  ///
  /// @param bValue is the binning value to be returned
  ///
  /// @return a one dimensional arrange
  const Range1D<ActsScalar>& range(BinningValue bValue) const;

  /// Return the N-dimension range
  const RangeXD<binValues, ActsScalar> range() const;

  /// Return an D-dimensional sub range according to the
  /// the given @param binValues
  template <unsigned int kSUBDIM>
  RangeXD<kSUBDIM, ActsScalar> range(
      const std::array<BinningValue, kSUBDIM>& binValues) const {
    RangeXD<kSUBDIM, ActsScalar> rRange;
    for (auto [i, v] : enumerate(binValues)) {
      rRange[i] = range(v);
    }
    return rRange;
  }

  /// Return the envelope
  const ExtentEnvelope& envelope() const;

  /// Access the minimum parameter
  ///
  /// @param bValue the binning identification
  ActsScalar min(BinningValue bValue) const { return m_range[bValue].min(); }

  /// Access the maximum parameter
  ///
  /// @param bValue the binning identification
  ActsScalar max(BinningValue bValue) const { return m_range[bValue].max(); }

  /// Access the maximum parameter
  ///
  /// @param bValue the binning identification
  ActsScalar medium(BinningValue bValue) const {
    return 0.5 * (m_range[bValue].min() + m_range[bValue].max());
  }

  /// Contains check
  ///
  /// @param rhs the extent that is check if it is contained
  /// @param bValue is the binning value, if set to binValues
  ///               the check on all is done
  ///
  /// @return true if the @param rhs is contained
  bool contains(const GeometricExtent& rhs,
                BinningValue bValue = binValues) const;

  /// Contains check
  ///
  /// @param rhs the extent that is check for intersection
  /// @param bValue is the binning value, if set to binValues
  ///               the check on all is done
  ///
  /// @return true if the @param rhs is contained
  bool intersects(const GeometricExtent& rhs,
                  BinningValue bVbVAlueal = binValues) const;

  /// Constraints check
  ///
  /// @param bValue is the binning value, if all the check on all is done
  bool constrains(BinningValue bValue = binValues) const;

  /// Convert to output stream for screen output
  /// @param sl [in,out] The output stream
  std::ostream& toStream(std::ostream& sl) const;

 private:
  std::bitset<binValues> m_constrains{0};
  RangeXD<binValues, ActsScalar> m_range;
  ExtentEnvelope m_envelope = zeroEnvelopes;
};

inline const RangeXD<binValues, ActsScalar> GeometricExtent::range() const {
  return m_range;
}

inline const ExtentEnvelope& GeometricExtent::envelope() const {
  return m_envelope;
}

/// Overload of << operator for std::ostream for debug output
std::ostream& operator<<(std::ostream& sl, const GeometricExtent& rhs);

}  // namespace Acts
