// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::HoughTransformSeeder::HoughTransformSeeder(
    ActsExamples::HoughTransformSeeder::Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("HoughTransformSeeder", lvl),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing space point input collections");
  }
  for (const auto& i : m_cfg.inputSpacePoints) {
    if (i.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collection");
  }

  if (m_cfg.gridConfig.rMax != m_cfg.seedFinderConfig.rMax) {
    throw std::invalid_argument("Inconsistent config rMax");
  }

  if (m_cfg.gridConfig.deltaRMax != m_cfg.seedFinderConfig.deltaRMax) {
    throw std::invalid_argument("Inconsistent config deltaRMax");
  }

  // Fill convenience variables
  m_step_x = (m_parMax[m_par_x] - m_parMin[m_par_x]) / m_imageSize_x;
  m_step_y = (m_parMax[m_par_y] - m_parMin[m_par_y]) / m_imageSize_y;
  for (unsigned i = 0; i <= m_imageSize_x; i++)
     m_bins_x.push_back(unquant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, i));
  for (unsigned i = 0; i <= m_imageSize_y; i++)
     m_bins_y.push_back(unquant(m_parMin[m_par_y], m_parMax[m_par_y], m_imageSize_y, i));
  
}

ActsExamples::ProcessCode ActsExamples::HoughTransformSeeder::execute(
   const AlgorithmContext& ctx) const {
   // construct the combined input container of space point pointers from all
   // configured input sources.
   // pre-compute the total size required so we only need to allocate once
   size_t nSpacePoints = 0;
   for (const auto& isp : m_cfg.inputSpacePoints) {
      nSpacePoints += ctx.eventStore.get<SimSpacePointContainer>(isp).size();
   }
   
   std::vector<const SimSpacePoint*> spacePointPtrs;
   spacePointPtrs.reserve(nSpacePoints);
   for (const auto& isp : m_cfg.inputSpacePoints) {
      for (const auto& spacePoint :
              ctx.eventStore.get<SimSpacePointContainer>(isp)) {
         // since the event store owns the space points, their pointers should be
         // stable and we do not need to create local copies.
         spacePointPtrs.push_back(&spacePoint);
      }
   }

  // run the seeding
  static thread_local SimSeedContainer seeds;
  seeds.clear();
  static thread_local decltype(finder)::State state;

  m_image = createImage(spacePointPtrs);
  
  for (unsigned y = 0; y < m_imageSize_y; y++)
     for (unsigned x = 0; x < m_imageSize_x; x++)
        if (passThreshold(m_image, x, y))
        {
           // Now find the set of hits for each such road JAAA TO DO need to know the layer!
           seeds.insert(a.end(), m_image(y, x).second); // JAA TO DO this is a collection of            
        }

  // extract proto tracks, i.e. groups of measurement indices, from tracks seeds
  size_t nSeeds = seeds.size();
  static thread_local ProtoTrackContainer protoTracks;
  protoTracks.clear();

  protoTracks.reserve(nSeeds);
  for (const auto& seed : seeds) {
    ProtoTrack& protoTrack = protoTracks.emplace_back();
    protoTrack.reserve(seed.sp().size());
    for (auto spacePointPtr : seed.sp()) {
      protoTrack.push_back(spacePointPtr->measurementIndex());
    }
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  ctx.eventStore.add(m_cfg.outputSeeds, SimSeedContainer{seeds});
  ctx.eventStore.add(m_cfg.outputProtoTracks, ProtoTrackContainer{protoTracks});
  return ActsExamples::ProcessCode::SUCCESS;
}



Image createLayerImage(unsigned layer, std::vector<const SimSpacePoint*> const & spacepoints) const {
{
    Image image(m_imageSize_y, m_imageSize_x);

    for (auto sp : spacepoints) 
    {
       if (sp->layer() != layer) continue; // JAAA

       // This scans over y (pT) because that is more efficient in memory, in C.
       // Unknown if firmware will want to scan over x instead.
       for (unsigned y_ = 0; y_ < m_imageSize_y; y_++)
       {
          unsigned y_bin_min = scale * y_;
          unsigned y_bin_max = scale * (y_ + 1);
          
          // Find the min/max x bins
          float r = TMath::Sqrt(sp.x()*sp.x()+sp.y()*sp.y());
          float phi = TMath::Atan2(sp.y()/sp.x());
          auto xBins = yToXBins(y_bin_min, y_bin_max, r, phi, sp.layer()); // JAAA
          
          // Update the image
          for (unsigned y = y_bin_min; y < y_bin_max; y++)
             for (unsigned x = xBins.first; x < xBins.second; x++)
             {
                image(y, x).first++;
                image(y, x).second.insert(hit);
             }
       }
    }
    
    return image;
}

HoughTransformSeeder::Image HoughTransformSeeder::createImage(std::vector<const SimSpacePoint*> const & spacepoints) const
{
    Image image(m_imageSize_y, m_imageSize_x);

    for (unsigned i = 0; i < m_nLayers; i++)
    {
        Image layerImage = createLayerImage(i, spacepoints);
        for (unsigned x = 0; x < m_imageSize_x; ++x)
            for (unsigned y = 0; y < m_imageSize_y; ++y)
                if (layerImage(y, x).first > 0)
                {
                    image(y, x).first++;
                    image(y, x).second.insert(layerImage(y, x).second.begin(), layerImage(y, x).second.end());
                }
    }

    return image;
}

bool HoughTransformSeeder::passThreshold(Image const & image, unsigned x, unsigned y) const
{
    // Pass window threshold
    if (x < 0 || image.size(1) < 0) return false;
    if (image(y, x.first < m_threshold) return false;

    // Pass local-maximum check
    if (m_localMaxWindowSize)
        for (int j = -m_localMaxWindowSize; j <= m_localMaxWindowSize; j++)
            for (int i = -m_localMaxWindowSize; i <= m_localMaxWindowSize; i++)
            {
                if (i == 0 && j == 0) continue;
                if (y + j < image.size(0) && x + i < image.size(1))
                {
                    if (image(y+j, x+i).first > image(y, x).first) return false;
                    if (image(y+j, x+i).first == image(y, x).first)
                    {
                        if (image(y+j, x+i).second.size() > image(y, x).second.size()) return false;
                        if (image(y+j, x+i).second.size() == image(y, x).second.size()
                            && j <= 0 && i <= 0) return false; // favor bottom-left (low phi, low neg q/pt)
                    }
                }
            }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Helpers


// Quantizes val, given a range [min, max) split into nSteps. Returns the bin below.
static inline int quant(double min, double max, unsigned nSteps, double val)
{
    return static_cast<int>((val - min) / (max - min) * nSteps);
}

// Returns the lower bound of the bin specified by step
static inline double unquant(double min, double max, unsigned nSteps, int step)
{
    return min + (max - min) * step / nSteps;
}

template <typename T>
static inline std::string to_string(std::vector<T> v)
{
    std::ostringstream oss;
    oss << "[";
    if (!v.empty())
    {
        std::copy(v.begin(), v.end()-1, std::ostream_iterator<T>(oss, ", "));
        oss << v.back();
    }
    oss << "]";
    return oss.str();
}


double HoughTransformSeeder::fieldCorrection(unsigned region, double qpt, double r)
{
   // can ultimately derive field corrections here!
   return 0;
}

double HoughTransformSeeder::yToX(double y, double r, double phi) const
{
   double d0 = 0; // d0 correction TO DO allow for this
   double x = asin(r * A * y - d0 / r) + phi;
    
   if (m_fieldCorrection) x += fieldCorrection(0, y, r);

   return x;
}

// Find the min/max x bins of the hit's line, in each y bin. Max is exclusive.
// Note this assumes yToX is monotonic. Returns {0, 0} if hit lies out of bounds.
std::pair<unsigned, unsigned> HoughTransformSeeder::yToXBins(size_t yBin_min, size_t yBin_max, double r, double phi, unsigned layer) const
{
    // Get float values
   double x_min = yToX(m_bins_y[yBin_min], r, phi);
   double x_max = yToX(m_bins_y[yBin_max], r, phi);
   if (x_min > x_max) std::swap(x_min, x_max);
   if (x_max < m_parMin[m_par_x] || x_min > m_parMax[m_par_x])
      return { 0, 0 }; // out of bounds
   
   // Get bins
   int x_bin_min = quant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, x_min);
   int x_bin_max = quant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, x_max) + 1; // exclusive

   // Extend bins
   unsigned extend = getExtension(yBin_min, layer);  // JAAAAAAAAAAAAAA TO DO LAYER
   x_bin_min -= extend;
   x_bin_max += extend;
   
   // Clamp bins
   if (x_bin_min < 0) x_bin_min = 0;
   if (x_bin_max > static_cast<int>(m_imageSize_x)) x_bin_max = m_imageSize_x;
   
   return { x_bin_min, x_bin_max };
}

// We allow variable extension based on the size of m_hitExtend_x. See comments below.
unsigned HoughTransformSeeder::getExtension(unsigned y, unsigned layer) const
{
   if (m_hitExtend_x.size() == m_nLayers) return m_hitExtend_x[layer];
   if (m_hitExtend_x.size() == m_nLayers * 2)
   {
      // different extension for low pt vs high pt, split in half but irrespective of sign
      // first nLayers entries of m_hitExtend_x is for low pt half, rest are for high pt half
      if (y < m_imageSize_y / 4 || y > 3 * m_imageSize_y / 4) return m_hitExtend_x[layer];
      return m_hitExtend_x[m_nLayers + layer];
   }
   return 0;
}

