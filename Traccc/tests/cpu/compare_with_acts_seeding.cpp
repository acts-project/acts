/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// io
#include "traccc/io/read_spacepoints.hpp"

// algorithms
#include "traccc/seeding/detail/seed_finding.hpp"
#include "traccc/seeding/detail/spacepoint_binning.hpp"

// tests
#include "tests/atlas_cuts.hpp"

// acts
#include <Acts/EventData/Seed.hpp>
#include <Acts/EventData/SpacePointContainer2.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/Seeding2/BroadTripletSeedFilter.hpp>
#include <Acts/Seeding2/CylindricalSpacePointGrid2.hpp>
#include <Acts/Seeding2/TripletSeeder.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/GridBinFinder.hpp>
#include <Acts/Utilities/Helpers.hpp>

// VecMem
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <algorithm>
#include <functional>
#include <limits>

class CompareWithActsSeedingTests
    : public ::testing::TestWithParam<
          std::tuple<std::string, std::string, unsigned int>> {};

// This defines the local frame test suite
TEST_P(CompareWithActsSeedingTests, Run) {

    // Skip this test, as the Acts and traccc seeders are not correctly
    // configured by it to behave the same way.
    // TODO: figure out what to do with this test.
    GTEST_SKIP() << "Acts and traccc seeders are out of sync";

    const std::string detector_file = std::get<0>(GetParam());
    const std::string hits_dir = std::get<1>(GetParam());
    const unsigned int event = std::get<2>(GetParam());

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Seeding Config
    traccc::seedfinder_config traccc_config;
    traccc_config.zMin = -1186.f * traccc::unit<float>::mm;
    traccc_config.zMax = 1186.f * traccc::unit<float>::mm;
    traccc_config.cotThetaMax = 7.40627f;
    traccc_config.deltaRMin = 1.f * traccc::unit<float>::mm;
    traccc_config.deltaRMax = 60.f * traccc::unit<float>::mm;
    traccc_config.sigmaScattering = 1.0f;
    traccc_config.deltaZMax = 1000000.f * traccc::unit<float>::mm;

    traccc::spacepoint_grid_config grid_config(traccc_config);

    traccc::seedfilter_config traccc_filter_config;

    // Declare algorithms
    traccc::host::details::spacepoint_binning sb{traccc_config, grid_config,
                                                 host_mr};
    traccc::host::details::seed_finding sf{traccc_config, traccc_filter_config,
                                           host_mr};

    // Read the hits from the relevant event file
    traccc::edm::spacepoint_collection::host spacepoints_per_event{host_mr};
    traccc::edm::measurement_collection::host measurements_per_event{host_mr};
    traccc::io::read_spacepoints(spacepoints_per_event, measurements_per_event,
                                 event, hits_dir);

    // Create Acts spacepoints from the traccc ones.
    Acts::SpacePointContainer2 actsSpacepoints{
        Acts::SpacePointColumns::PackedXY | Acts::SpacePointColumns::PackedZR |
        Acts::SpacePointColumns::VarianceZ |
        Acts::SpacePointColumns::VarianceR |
        Acts::SpacePointColumns::CopyFromIndex};
    actsSpacepoints.reserve(
        static_cast<unsigned int>(spacepoints_per_event.size()));
    for (size_t i = 0; i < spacepoints_per_event.size(); ++i) {
        // A convenience proxy for the traccc spacepoint:
        traccc::edm::spacepoint tsp = spacepoints_per_event[i];
        // Create the corresponding Acts spacepoint.
        auto asp = actsSpacepoints.createSpacePoint();
        asp.xy() = {static_cast<float>(tsp.x()), static_cast<float>(tsp.y())};
        asp.zr() = {static_cast<float>(tsp.z()),
                    static_cast<float>(tsp.radius())};
        asp.varianceR() = static_cast<float>(tsp.radius_variance());
        asp.varianceZ() = static_cast<float>(tsp.z_variance());
        asp.copyFromIndex() = static_cast<unsigned int>(i);
    }

    /*--------------------------------
      TRACCC seeding
      --------------------------------*/

    auto internal_spacepoints_per_event =
        sb(vecmem::get_data(spacepoints_per_event));
    auto traccc_seeds = sf(vecmem::get_data(spacepoints_per_event),
                           internal_spacepoints_per_event);

    /*--------------------------------
      ACTS seeding
      --------------------------------*/

    // Create the grid. This is using a CylindricalSpacePointGrid which is
    // a 3D grid (phi, z, radius). We need to pass a config and an option
    // objects
    Acts::CylindricalSpacePointGrid2::Config gridConf;
    gridConf.minPt = traccc_config.minPt;
    gridConf.rMax = traccc_config.rMax;
    gridConf.rMin = 0;
    gridConf.zMax = traccc_config.zMax;
    gridConf.zMin = traccc_config.zMin;
    gridConf.deltaRMax = traccc_config.deltaRMax;
    gridConf.cotThetaMax = traccc_config.cotThetaMax;
    gridConf.impactMax = traccc_config.impactMax;
    gridConf.phiMin = traccc_config.phiMin;
    gridConf.phiMax = traccc_config.phiMax;
    // The Bin Ednges allow the user to define a variable binnin in the z or
    // radius azis. If no value is provided, the values of zMin/zMax and
    // rMin/rMax will be used For the phi axis the number of binnings are
    // compute inside the createGrid function
    gridConf.zBinEdges = {};
    gridConf.rBinEdges = {};
    gridConf.phiBinDeflectionCoverage = traccc_config.phiBinDeflectionCoverage;
    gridConf.bFieldInZ = traccc_config.bFieldInZ;
    gridConf.bottomBinFinder.emplace(1, std::vector<std::pair<int, int>>{}, 0);
    gridConf.topBinFinder.emplace(1, std::vector<std::pair<int, int>>{}, 0);
    gridConf.navigation[0ul] = {};
    gridConf.navigation[1ul] = {};
    gridConf.navigation[2ul] = {};

    Acts::DoubletSeedFinder::Config bottomDoubletFinderConfig;
    bottomDoubletFinderConfig.candidateDirection = Acts::Direction::Backward();
    bottomDoubletFinderConfig.deltaRMin = traccc_config.deltaRMin;
    bottomDoubletFinderConfig.deltaRMax = traccc_config.deltaRMax;
    bottomDoubletFinderConfig.deltaZMin = 0.f;
    bottomDoubletFinderConfig.deltaZMax = traccc_config.deltaZMax;
    bottomDoubletFinderConfig.impactMax = traccc_config.impactMax;
    bottomDoubletFinderConfig.collisionRegionMin =
        traccc_config.collisionRegionMin;
    bottomDoubletFinderConfig.collisionRegionMax =
        traccc_config.collisionRegionMax;
    bottomDoubletFinderConfig.cotThetaMax = traccc_config.cotThetaMax;
    bottomDoubletFinderConfig.minPt = traccc_config.minPt;
    auto bottomDoubletFinder =
        Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
            bottomDoubletFinderConfig, traccc_config.bFieldInZ));

    Acts::DoubletSeedFinder::Config topDoubletFinderConfig =
        bottomDoubletFinderConfig;
    topDoubletFinderConfig.candidateDirection = Acts::Direction::Forward();
    auto topDoubletFinder =
        Acts::DoubletSeedFinder::create(Acts::DoubletSeedFinder::DerivedConfig(
            topDoubletFinderConfig, traccc_config.bFieldInZ));

    Acts::TripletSeedFinder::Config tripletFinderConfig;
    tripletFinderConfig.useStripInfo = false;
    tripletFinderConfig.sortedByCotTheta = true;
    tripletFinderConfig.minPt = traccc_config.minPt;
    tripletFinderConfig.sigmaScattering = traccc_config.sigmaScattering;
    tripletFinderConfig.radLengthPerSeed = traccc_config.radLengthPerSeed;
    tripletFinderConfig.impactMax = traccc_config.impactMax;
    auto tripletFinder =
        Acts::TripletSeedFinder::create(Acts::TripletSeedFinder::DerivedConfig(
            tripletFinderConfig, traccc_config.bFieldInZ));

    Acts::BroadTripletSeedFilter::Config filterConfig;
    filterConfig.deltaInvHelixDiameter =
        traccc_filter_config.deltaInvHelixDiameter;
    filterConfig.deltaRMin = traccc_filter_config.deltaRMin;
    filterConfig.compatSeedWeight = traccc_filter_config.compatSeedWeight;
    filterConfig.impactWeightFactor = traccc_filter_config.impactWeightFactor;
    filterConfig.maxSeedsPerSpM = traccc_config.maxSeedsPerSpM;
    filterConfig.compatSeedLimit = traccc_filter_config.compatSeedLimit;

    // We fill the grid with the space points
    Acts::CylindricalSpacePointGrid2 grid(gridConf);
    for (unsigned int i = 0; i < actsSpacepoints.size(); ++i) {
        const auto& sp = actsSpacepoints[i];
        float phi = std::atan2(sp.xy()[1], sp.xy()[0]);
        grid.insert(i, phi, sp.zr()[0], sp.zr()[1]);
    }

    // Get the traccc axes:
    //  0 -> phi
    //  1 -> z
    traccc::axis2::circular axis0 = internal_spacepoints_per_event.axis_p0();
    traccc::axis2::regular axis1 = internal_spacepoints_per_event.axis_p1();

    const auto& gridAxes = grid.grid().axes();
    auto acts_axis0 = gridAxes[0];
    auto acts_axis1 = gridAxes[1];

    EXPECT_EQ(axis0.bins(), acts_axis0->getNBins());
    EXPECT_EQ(axis1.bins(), acts_axis1->getNBins());

    auto axis0_borders = axis0.all_borders();
    auto axis1_borders = axis1.all_borders();
    auto acts_axis0_borders = acts_axis0->getBinEdges();
    auto acts_axis1_borders = acts_axis1->getBinEdges();

    EXPECT_EQ(axis0_borders.size(), axis0.bins() + 1);
    EXPECT_EQ(axis1_borders.size(), axis1.bins() + 1);

    EXPECT_EQ(axis0_borders.size(), acts_axis0_borders.size());
    EXPECT_EQ(axis1_borders.size(), acts_axis1_borders.size());

    for (unsigned int i = 0; i < axis0.bins() + 1; ++i) {
        EXPECT_NEAR(axis0_borders[i], acts_axis0_borders[i], 0.01);
    }
    for (unsigned int i = 0; i < axis1.bins() + 1; ++i) {
        EXPECT_NEAR(axis1_borders[i], acts_axis1_borders[i], 0.01);
    }

    // Compute radius Range for the middle space point candidate
    // we rely on the fact the grid is storing the proxies
    // with a sorting in the radius
    float minRange = std::numeric_limits<float>::max();
    float maxRange = std::numeric_limits<float>::lowest();
    for (const Acts::CylindricalSpacePointGrid2::BinType& coll : grid.grid()) {
        if (coll.empty()) {
            continue;
        }
        const auto firstEl = actsSpacepoints.at(coll.front());
        const auto lastEl = actsSpacepoints.at(coll.back());
        minRange = std::min(firstEl.zr()[1], minRange);
        maxRange = std::max(lastEl.zr()[1], maxRange);
    }

    const Acts::Range1D<float> rMiddleSPRange(
        std::floor(minRange / 2) * 2 +
            10.f * static_cast<float>(Acts::UnitConstants::mm),
        std::floor(maxRange / 2) * 2 -
            10.f * static_cast<float>(Acts::UnitConstants::mm));

    // We define the state and the seed container
    Acts::SeedContainer2 actsSeeds;
    actsSeeds.assignSpacePointContainer(actsSpacepoints);

    // Run the Acts seeding
    Acts::BroadTripletSeedFilter::State filterState;
    Acts::BroadTripletSeedFilter::Cache filterCache;
    Acts::BroadTripletSeedFilter seedFilter(
        filterConfig, filterState, filterCache, Acts::getDummyLogger());
    Acts::TripletSeeder actsSeedfinder;
    Acts::TripletSeeder::Cache cache;
    for (const auto [bottom, middle, top] : grid.binnedGroup()) {

        // Collect all the bottom and top spacepoints from the grid bins
        // provided by "bottom" and "top".
        std::vector<unsigned int> bottomIndices;
        for (std::size_t binIndex : bottom) {
            const Acts::CylindricalSpacePointGrid2::BinType& bottomBin =
                grid.at(binIndex);
            bottomIndices.insert(bottomIndices.end(), bottomBin.begin(),
                                 bottomBin.end());
        }
        std::vector<unsigned int> topIndices;
        for (std::size_t binIndex : top) {
            const Acts::CylindricalSpacePointGrid2::BinType& topBin =
                grid.at(binIndex);
            topIndices.insert(topIndices.end(), topBin.begin(), topBin.end());
        }

        // Loop over all the middle spacepoints in the bin provided by "middle".
        for (unsigned int middleIndex : grid.at(middle)) {

            // Only consider middle spacepoints in the defined radius range
            const auto middleSP = actsSpacepoints[middleIndex].asConst();
            if (middleSP.zr()[1] < rMiddleSPRange.min() ||
                middleSP.zr()[1] > rMiddleSPRange.max()) {
                continue;
            }

            // Find the seeds for this specific middle spacepoint.
            auto bottomSpacepoints =
                actsSpacepoints.subset(bottomIndices).asConst();
            auto topSpacepoints = actsSpacepoints.subset(topIndices).asConst();
            actsSeedfinder.createSeedsFromGroup(
                cache, *bottomDoubletFinder, *topDoubletFinder, *tripletFinder,
                seedFilter, actsSpacepoints, bottomSpacepoints, middleSP,
                topSpacepoints, actsSeeds);
        }
    }

    // We have created seeds of proxies to space point at this point
    // From a proxy we can retrieve the original space point object with
    // the externalSpacePoint() method

    // Count the number of matching seeds
    std::size_t n_matched_acts_seeds = 0u;
    for (traccc::edm::seed_collection::host::size_type i = 0u;
         i < traccc_seeds.size(); ++i) {
        const auto traccc_seed = traccc_seeds.at(i);
        // Try to find the same Acts seed.
        auto it = std::find_if(
            actsSeeds.begin(), actsSeeds.end(), [&](const auto& acts_seed) {
                return ((traccc_seed.bottom_index() ==
                         actsSpacepoints[acts_seed.spacePointIndices()[0]]
                             .copyFromIndex()) &&
                        (traccc_seed.middle_index() ==
                         actsSpacepoints[acts_seed.spacePointIndices()[1]]
                             .copyFromIndex()) &&
                        (traccc_seed.top_index() ==
                         actsSpacepoints[acts_seed.spacePointIndices()[2]]
                             .copyFromIndex()));
            });

        if (it != actsSeeds.end()) {
            ++n_matched_acts_seeds;
        }
    }

    // Ensure that ACTS and traccc give the same result
    // @TODO Uncomment the line below once acts-project/acts#2132 is merged
    // EXPECT_EQ(seeds.size(), seedVector.size());
    EXPECT_NEAR(static_cast<double>(traccc_seeds.size()),
                static_cast<double>(n_matched_acts_seeds),
                static_cast<double>(traccc_seeds.size()) * 0.004);
    EXPECT_GT(static_cast<double>(n_matched_acts_seeds) /
                  static_cast<double>(traccc_seeds.size()),
              0.995);
}

INSTANTIATE_TEST_SUITE_P(
    SeedingValidation, CompareWithActsSeedingTests,
    ::testing::Values(std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 0),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 1),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 2),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 3),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 4),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 5),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 6),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 7),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 8),
                      std::make_tuple("tml_detector/trackml-detector.csv",
                                      "tml_full/ttbar_mu200/", 9)));
