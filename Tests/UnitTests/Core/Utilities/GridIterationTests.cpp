// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridIterator.hpp"

#include <array>

namespace Acts::Test {

  BOOST_AUTO_TEST_CASE(grid_iteration_test_1d_global) {
    const std::size_t nBins = 10ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::Grid<double, Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis)) );

    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins);
    BOOST_CHECK_EQUAL(grid.size(true), nBins + 2ul);

    const std::array<std::size_t, 1ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);

    Acts::GridGlobalIterator<double, Acts::detail::EquidistantAxis> gridStart = grid.begin();
    Acts::GridGlobalIterator<double, Acts::detail::EquidistantAxis> gridStop = grid.end();
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(true));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_2d_global) {
    const std::size_t nBins = 10ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::Grid<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis)) );

    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul));

    const std::array<std::size_t, 2ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);

    Acts::GridGlobalIterator<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> gridStart = grid.begin();
    Acts::GridGlobalIterator<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> gridStop = grid.end();
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(true));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_3d_global) {
    const std::size_t nBins = 10ul;
    const std::size_t nBinsZ = 20ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis zAxis(0, 100, nBinsZ);
    Acts::Grid<double,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)) );

    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins * nBinsZ);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul) * (nBinsZ + 2ul));

    const std::array<std::size_t, 3ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[2ul], nBinsZ);

    Acts::GridGlobalIterator<double,
			     Acts::detail::EquidistantAxis,
			     Acts::detail::EquidistantAxis,
			     Acts::detail::EquidistantAxis> gridStart = grid.begin();
    Acts::GridGlobalIterator<double,
			     Acts::detail::EquidistantAxis,
			     Acts::detail::EquidistantAxis,
			     Acts::detail::EquidistantAxis> gridStop = grid.end();
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(true));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_1d_local) {
    const std::size_t nBins = 10ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::Grid<double, Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis)) );

    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins);
    BOOST_CHECK_EQUAL(grid.size(true), nBins + 2ul);

    const std::array<std::size_t, 1ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);

    std::array<std::vector<std::size_t>, 1ul> navigation;
    navigation[0ul].resize(nBins);
    std::iota(navigation[0ul].begin(), navigation[0ul].end(), 1);
    
    Acts::GridLocalIterator<double, Acts::detail::EquidistantAxis> gridStart = grid.begin(navigation);
    Acts::GridLocalIterator<double, Acts::detail::EquidistantAxis> gridStop = grid.end(navigation);
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(false));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_2d_local) {
    const std::size_t nBins = 10ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::Grid<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis)) );
    
    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul));

    const std::array<std::size_t, 2ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);

    std::array<std::vector<std::size_t>, 2ul> navigation;
    navigation[0ul].resize(nBins);
    navigation[1ul].resize(nBins);
    for (std::size_t i(0ul); i<2ul; ++i) {
      std::iota(navigation[i].begin(), navigation[i].end(), 1);
    }
    
    Acts::GridLocalIterator<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> gridStart = grid.begin(navigation);
    Acts::GridLocalIterator<double, Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis> gridStop = grid.end(navigation);
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(false));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_3d_local) {
    const std::size_t nBins = 10ul;
    const std::size_t nBinsZ = 20ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis zAxis(0, 100, nBinsZ);
    Acts::Grid<double,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)) );
    
    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins * nBinsZ);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul) * (nBinsZ + 2ul));

    const std::array<std::size_t, 3ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[2ul], nBinsZ);

    std::array<std::vector<std::size_t>, 3ul> navigation;
    navigation[0ul].resize(nBins);
    navigation[1ul].resize(nBins);
    navigation[2ul].resize(nBinsZ);
    for (std::size_t i(0ul); i<3ul; ++i) {
      std::iota(navigation[i].begin(), navigation[i].end(), 1);
    }
    
    Acts::GridLocalIterator<double,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis> gridStart = grid.begin(navigation);
    Acts::GridLocalIterator<double,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis> gridStop = grid.end(navigation);
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(false));
  }
  
  BOOST_AUTO_TEST_CASE(grid_iteration_test_3d_local_custom_navigation) {
    const std::size_t nBins = 10ul;
    const std::size_t nBinsZ = 20ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis zAxis(0, 100, nBinsZ);
    Acts::Grid<double,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis)) );
    
    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins * nBinsZ);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul) * (nBinsZ + 2ul));

    const std::array<std::size_t, 3ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[2ul], nBinsZ);

    std::array<std::vector<std::size_t>, 3ul> navigation;
    navigation[0ul] = {1ul, 5ul, 3ul, 2ul, 9ul, 10ul, 4ul, 6ul, 8ul, 7ul};
    navigation[1ul] = {6ul, 8ul, 7ul, 1ul, 5ul, 3ul, 2ul, 9ul, 10ul, 4ul};
    navigation[2ul] = {1ul, 5ul, 3ul, 2ul, 9ul, 10ul, 4ul, 6ul, 8ul, 7ul, 11ul, 15ul, 13ul, 12ul, 19ul, 20ul, 14ul, 16ul, 18ul, 17ul};
    
    Acts::GridLocalIterator<double,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis> gridStart = grid.begin(navigation);
    Acts::GridLocalIterator<double,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis> gridStop = grid.end(navigation);
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }
    BOOST_CHECK_EQUAL(numIterations, grid.size(false));
  }

  BOOST_AUTO_TEST_CASE(grid_iteration_test_5d_local_custom_subnavigation) {
    const std::size_t nBins = 10ul;
    const std::size_t nBinsZ = 20ul;
    const std::size_t nBinsJK = 5ul;
    Acts::detail::EquidistantAxis xAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis yAxis(0, 100, nBins);
    Acts::detail::EquidistantAxis zAxis(0, 100, nBinsZ);
    Acts::detail::EquidistantAxis jAxis(0, 100, nBinsJK);
    Acts::detail::EquidistantAxis kAxis(0, 100, nBinsJK);
    Acts::Grid<double,
               Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis,
	       Acts::detail::EquidistantAxis,
               Acts::detail::EquidistantAxis,
               Acts::detail::EquidistantAxis> grid( std::make_tuple(std::move(xAxis), std::move(yAxis), std::move(zAxis), std::move(jAxis), std::move(kAxis)) );


    // test general properties
    BOOST_CHECK_EQUAL(grid.size(false), nBins * nBins * nBinsZ * nBinsJK * nBinsJK);
    BOOST_CHECK_EQUAL(grid.size(true), (nBins + 2ul) * (nBins + 2ul) * (nBinsZ + 2ul) * (nBinsJK + 2ul) * (nBinsJK + 2ul));
    
    const std::array<std::size_t, 5ul> numLocalBins = grid.numLocalBins();
    BOOST_CHECK_EQUAL(numLocalBins[0ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[1ul], nBins);
    BOOST_CHECK_EQUAL(numLocalBins[2ul], nBinsZ);
    BOOST_CHECK_EQUAL(numLocalBins[3ul], nBinsJK);
    BOOST_CHECK_EQUAL(numLocalBins[4ul], nBinsJK);

    // Iterate only on a few bins
    std::array<std::vector<std::size_t>, 5ul> navigation;
    navigation[0ul] = {1ul, 5ul, 3ul, 2ul, 9ul, 10ul, 4ul, 6ul, 8ul, 7ul};
    navigation[1ul] = {6ul, 8ul, 7ul, 1ul};
    navigation[2ul] = {1ul, 5ul};
    navigation[3ul] = {5ul, 3ul, 2ul};
    navigation[4ul] = {2ul};

    Acts::GridLocalIterator<double,
                            Acts::detail::EquidistantAxis,
                            Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
                            Acts::detail::EquidistantAxis> gridStart = grid.begin(navigation);
    Acts::GridLocalIterator<double,
                            Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
			    Acts::detail::EquidistantAxis,
                            Acts::detail::EquidistantAxis,
                            Acts::detail::EquidistantAxis> gridStop = grid.end(navigation);
    std::size_t numIterations = 0ul;
    for (; gridStart != gridStop; ++gridStart) {
      ++numIterations;
    }

    std::size_t expectedIterations = 1ul;
    for (std::size_t i(0ul); i<5ul; ++i) {
      expectedIterations *= navigation[i].size();
    }      

    BOOST_CHECK_EQUAL(numIterations, expectedIterations);
  }
  
}
