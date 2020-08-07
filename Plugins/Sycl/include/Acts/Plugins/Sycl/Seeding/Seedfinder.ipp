#include <cmath>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <exception>
#include <boost/range/adaptors.hpp>
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.hpp"

namespace Acts::Sycl {
template <typename external_spacepoint_t>
Seedfinder<external_spacepoint_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(std::move(config)) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;
  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);

  // catch asynchronous exceptions
  auto exception_handler = [] (cl::sycl::exception_list exceptions) {
  for (std::exception_ptr const& e : exceptions) {
      try {
        std::rethrow_exception(e);
      } catch(cl::sycl::exception const& e) {
        std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
      }
    }
  };

  // create queue with costum device selector
  m_queue = cl::sycl::queue(nvidia_selector(), exception_handler);
}

template <typename external_spacepoint_t>
template <typename sp_range_t>
std::vector<Acts::Seed<external_spacepoint_t>>
Seedfinder<external_spacepoint_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  std::vector<float> offloadBottomSPs;
  std::vector<float> offloadMiddleSPs;
  std::vector<float> offloadTopSPs;

  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> bottomSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> middleSPvec;
  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> topSPvec;

  for(auto SP: bottomSPs) {
    bottomSPvec.push_back(SP);
    offloadBottomSPs.insert(offloadBottomSPs.end(),
                            {SP->x(), SP->y(), SP->z(), SP->radius(),
                             SP->varianceR(), SP->varianceZ()});
  }

  for(auto SP: middleSPs) {
    middleSPvec.push_back(SP);
    offloadMiddleSPs.insert(offloadMiddleSPs.end(),
                            {SP->x(), SP->y(), SP->z(), SP->radius(),
                             SP->varianceR(), SP->varianceZ()});
  }

  for(auto SP: topSPs) {
    topSPvec.push_back(SP);
    offloadTopSPs.insert(offloadTopSPs.end(),
                         {SP->x(), SP->y(), SP->z(), SP->radius(),
                          SP->varianceR(), SP->varianceZ()});
  }

  // sort top space points BY RADIUS for later filter algorithm

  // std::sort(topSPvec.begin(), topSPvec.end(), [](auto sp1, auto sp2){return sp1->radius() < sp2->radius();});
  // offloadTopSPs.reserve(topSPvec.size() * int(eSP));
  // for(auto SP: topSPvec) {
  //   offloadTopSPs.insert(offloadTopSPs.end(),
  //                        {SP->x(), SP->y(), SP->z(), SP->radius(),
  //                         SP->varianceR(), SP->varianceZ()});
  // }

  // turns out they are already sorted

  const int numBottomSPs = bottomSPvec.size();
  const int numMiddleSPs = middleSPvec.size();
  const int numTopSPs = topSPvec.size();

  // order is important because of enum order (see eConfigData)
  std::vector<float> offloadConfigData = {
    m_config.deltaRMin,
    m_config.deltaRMax,
    m_config.cotThetaMax,
    m_config.collisionRegionMin,
    m_config.collisionRegionMax,
    m_config.maxScatteringAngle2,
    m_config.sigmaScattering,
    m_config.minHelixDiameter2,
    m_config.pT2perRadius,
    m_config.seedFilter->getSeedFilterConfig().deltaInvHelixDiameter,
    m_config.seedFilter->getSeedFilterConfig().impactWeightFactor,
    m_config.seedFilter->getSeedFilterConfig().deltaRMin,
    m_config.seedFilter->getSeedFilterConfig().compatSeedWeight,
    float(m_config.seedFilter->getSeedFilterConfig().compatSeedLimit),
    m_config.impactMax
  };

  std::vector<std::vector<int>> seedIndices;
  std::vector<std::vector<float>>  seedWeight;

  offloadComputations(m_queue,
                      offloadConfigData,
                      offloadBottomSPs,
                      offloadMiddleSPs,
                      offloadTopSPs,
                      seedIndices,
                      seedWeight
  );


  for(int mi = 0; mi < numMiddleSPs && mi < seedIndices.size() && mi < seedWeight.size(); ++mi) {
    std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
        seedsPerSPM;
    for(int j = 0; j*2+1 < seedIndices[mi].size() && 
                    seedIndices[mi][2*j] != -1 &&
                    seedIndices[mi][2*j+1] != -1 ; ++j){
      auto& bottomSP = *(bottomSPvec[seedIndices[mi][2*j]]);
      auto& middleSP = *(middleSPvec[mi]);
      auto& topSP =    *(topSPvec[seedIndices[mi][2*j+1]]);
      float weight =   seedWeight[mi][j];

      seedsPerSPM.emplace_back(std::make_pair(weight, std::make_unique<const InternalSeed<external_spacepoint_t>>(
                        bottomSP, middleSP, topSP, 0)));

    }
     m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSPM, outputVec);
  }
  return outputVec;
}
}  // namespace Acts::Sycl