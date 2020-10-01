// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

namespace ActsFatras {

BOOST_AUTO_TEST_SUITE(Digitization)

BOOST_AUTO_TEST_CASE(CellChannelCluster) {
  // A single cell - loc0
  Cell cell0(5, 10.5);
  // One-dimensional channel
  Channel<Acts::eBoundLoc0> ch0 = {{cell0}, 0.5};

  // A second cell - loc1
  Cell cell1(6, -0.45);
  // Two-dimensional channel
  Channel<Acts::eBoundLoc0, Acts::eBoundLoc1> ch01 = {{cell0, cell1}, 0.75};

  // Creating a cluster
  Acts::ParameterSet<Acts::BoundIndices, Acts::eBoundLoc0,
                     Acts::eBoundLoc1>::ParametersVector pVec;
  pVec << cell0.second, cell1.second;
  Acts::ParameterSet<Acts::BoundIndices, Acts::eBoundLoc0,
                     Acts::eBoundLoc1>::CovarianceMatrix pCov;
  pCov << 0.1, 0., 0.3, 0.;
  Acts::ParameterSet<Acts::BoundIndices, Acts::eBoundLoc0, Acts::eBoundLoc1>
      pSet(std::move(pCov), pVec);

  Cluster<Acts::eBoundLoc0, Acts::eBoundLoc1> cluster(std::move(pSet), {1, 1},
                                                      {ch01});
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsFatras
