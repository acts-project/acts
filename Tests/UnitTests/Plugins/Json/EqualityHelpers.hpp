// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningData.hpp"

#include <iostream>

namespace Acts {

/// Check whether the BinningData objects are equal
///
/// @param ba The first BinningData object
/// @param bb The second BinningData object
/// @param tolerance a tolerance parameter
///
/// @return a boolean
inline static bool isEqual(const BinningData& ba, const BinningData& bb,
                           float tolerance) {
  bool equalBool = (ba.type == bb.type) and (ba.option == bb.option) and
                   (ba.binvalue == bb.binvalue) and (ba.zdim == bb.zdim) and
                   (ba.subBinningAdditive == bb.subBinningAdditive);

  BOOST_CHECK(equalBool);
  bool equalRange = (std::abs(ba.min - bb.min) < tolerance) and
                    (std::abs(ba.max - bb.max) < tolerance) and
                    (std::abs(ba.step - bb.step) < tolerance);

  BOOST_CHECK(equalRange);
  bool euqalStructure =
      (ba.subBinningData != nullptr)
          ? isEqual(*ba.subBinningData, *bb.subBinningData, tolerance)
          : (bb.subBinningData == nullptr);

  BOOST_CHECK(euqalStructure);

  bool equalBoundaries = (ba.boundaries().size() == bb.boundaries().size());
  if (equalBoundaries) {
    for (size_t ib = 0; ib < ba.boundaries().size(); ++ib) {
      equalBoundaries =
          (std::abs(ba.boundaries()[ib] - bb.boundaries()[ib]) < tolerance);
      if (not equalBoundaries) {
        break;
      }
    }
  }
  BOOST_CHECK(equalBoundaries);

  return equalBool and equalRange and euqalStructure;
}

/// Check whether the BinUtility objects are equal
///
/// @param ba The first BinUtility object
/// @param bb the second BinUtility object
/// @param tolerance a tolerance parameter
///
/// @return a bollean if equal
inline static bool isEqual(const BinUtility& ba, const BinUtility& bb,
                           float tolerance) {
  bool equal = (ba.binningData().size() == bb.binningData().size());
  BOOST_CHECK(equal);
  if (equal) {
    for (size_t ib = 0; ib < ba.binningData().size(); ++ib) {
      equal = isEqual(ba.binningData()[ib], bb.binningData()[ib], tolerance);
      BOOST_CHECK(equal);
    }
  }
  return equal;
}

/// check whether Extnet objects are equal - with tolerance
///
/// @param ea the first extent object
/// @param eb the second extent object
/// @param tolerance the tolerance parameter
///
/// @return bool for euqal
inline static bool isEqual(const Acts::Extent& ea, const Acts::Extent& eb,
                           Acts::ActsScalar tolerance = 0.) {
  bool equalConstrains = true;
  bool equalRange = true;
  for (auto& bVal : s_binningValues) {
    equalConstrains =
        equalConstrains and (ea.constrains(bVal) == eb.constrains(bVal));
    BOOST_CHECK(equalConstrains);
    if (ea.constrains(bVal) and eb.constrains(bVal)) {
      equalRange =
          equalRange and std::abs(ea.min(bVal) - eb.min(bVal)) < tolerance;
      equalRange =
          equalRange and std::abs(ea.max(bVal) - eb.max(bVal)) < tolerance;
      BOOST_CHECK(equalRange);
    }
  }
  BOOST_CHECK(equalConstrains);
  BOOST_CHECK(equalRange);
  return equalRange and equalConstrains;
}

}  // namespace Acts