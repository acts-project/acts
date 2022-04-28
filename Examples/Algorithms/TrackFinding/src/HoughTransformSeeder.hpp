// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// 
// @file HoughTransformSeeder.hpp
// @author Riley Xu then modified to ACTS by Jahred Adelman
// @brief Implements track-seeding using a Hough transform.
// 
// Using the Lorentz force equation, one can relate the phi of a track and the
// coordinate of a single hit:
// 
//      A * q / pT = sin(phi_track - phi_hit) / r
// 
// where
//      A   : 3 * 10^-4 GeV / (c*mm*e) for the ATLAS B field (to be configured)
//      q   : charge of the particle
//      pT  : transverse momentum
//      r   : cylindrical radius of the hit from the beamline
//      phi : in radians
// 
// Here, q/pT and phi_track are unknown. This equation forms a line in q/pT vs
// phi_track space. Each hit will have its own line based on its phi and r.
// However, note that hits belonging to the same track will have lines that
// intersect at the track's q/pT and phi. In this manner, we can conduct pattern
// -matching by looking for intersections of these pT-phi lines.
// 
// To easily find intersections, we first pixelate (equivalently, we make a 2d
// histogram from) the graph of all the hit's lines in q/pT vs phi_track space.
// We then apply a convolution (i.e. a scanning window) to pick out points with
// multiple lines going through them. These points become our roads.
// 
// In principle the Hough transform can be used for an entire region (i.e. .2 phi x .2 eta) or larger.
// However this can lead to an excessive number of hits/lines in the transform image, leading
// to spurious intersections. Instead, we can use multiple transforms that each cover a slice
// in z0, and simply combine all the roads found.
// 
// References:
//      Martensson Thesis: http://uu.diva-portal.org/smash/get/diva2:1341509/FULLTEXT01.pdf
// 

// We adopt the following nomenclature within this class:
//      image: The 'graph' in q/pT vs phi_track space, filled with a line calculated as above for each hit.
//      point: A specific q/pT and phi_track bin in the above image; i.e. what is normally called a pixel
//             but I don't want to confuse this with the detector type. A point's value is the number of
//             lines that go through it.
//
// For the first iteration, x refers to phi_track, and y refers to q/pT, although
// this should remain flexible. These are set via the variables m_par_x and m_par_y.
//
// NOTE: We store the image in graph sense and not computer-science sense. That is,
// the row-index is y. The y-axis still points downwards, so that y=0 represents the
// lowest bin.
//      image[y=0][x=0]      : lowest q/pT and lowest phi_track bin
//      image[y=size-1][x=0] : highest q/pT and lowest phi_track bin
//


#pragma once


#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"

#include <string>
#include <vector>
#include <utility>
#include <unordered_set>



namespace ActsExamples {

/// Construct track seeds from space points.
class HoughTransformSeeder final : public BareAlgorithm {
 public:
  struct Config {
    /// Input space point collections.
    ///
    /// We allow multiple space point collections to allow different parts of
    /// the detector to use different algorithms for space point construction,
    /// e.g. single-hit space points for pixel-like detectors or double-hit
    /// space points for strip-like detectors.
    std::vector<std::string> inputSpacePoints;
    /// Output track seed collection.
    std::string outputSeeds;
    /// Output proto track collection.
    std::string outputProtoTracks;

    Acts::SeedFilterConfig seedFilterConfig;
    Acts::SeedfinderConfig<SimSpacePoint> seedFinderConfig;
    Acts::SpacePointGridConfig gridConfig;
  };

  /// Construct the seeding algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  HoughTransformSeeder(Config cfg, Acts::Logging::Level lvl);

  /// Run the seeding algorithm.
  ///
  /// @param txt is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

  /// Const access to the config
  const Config& config() const { return m_cfg; }

  double getMinX() const { return m_parMin[m_par_x]; }
  double getMaxX() const { return m_parMax[m_par_x]; }
  double getMinY() const { return m_parMin[m_par_y]; }
  double getMaxY() const { return m_parMax[m_par_y]; }
  unsigned getThreshold() const { return m_threshold[m_threshold.size() / 2]; }
  int getSubRegion() const { return m_subRegion; }

  static constexpr double kA = 0.0003; // Assume B = 2T constant. TODO Make this a parameter
 
  // Apply correction due to B != 2T everywhere. This correction should be ADDED to
  // phi_track. // for now this does nothing!
  static double fieldCorrection(unsigned region, double y, double r); 
  
  double yToX(double y, double r, double phi) const;

  typedef Vector2D<std::pair<int, std::unordered_set<const SimSpacePoint*>>> Image;
  // An image is a 2d array of points, where each point has a value.
  // The value starts as the number of hit layers, but can change with effects
  // like a convolution. Also stored are all hits that contributed to each bin.
  // Size m_imageSize_y * m_imageSize_x. (NOTE y is row coordinate)
  Image const & getImage() const { return m_image; } // Returns the image generated from the last call of getSeeds()
  
  

 private:
  Config m_cfg;
  int m_subRegion = 0; // -1 for entire region (no slicing), slicing is to be added TODO

  // Trace each hit that goes in a bin
  // Disabling this will save memory/time since each bin doesn't have to store all its hits
  // but the roads created won't have hits from convolution, etc.
  bool m_traceHits = true;

  bool m_fieldCorrection = true;

  // === Image ===
  float m_xMin = 0.3; // minphi
  float m_xMax = 0.5; // maxphi
  float m_yMin = -1.0; // min q/pt, -1/1 GeV
  float m_yMax = 1.0; // max q/pt, +1/1 GeV

  unsigned m_imageSize_x = 0; // i.e. number of bins in phi_track
  unsigned m_imageSize_y = 0; // i.e. number of bins in q/pT
  
  std::vector<unsigned> m_hitExtend_x; // Hit lines will fill extra bins in x by this amount on each side, size == nLayers

  // === Seeds ===
  int m_threshold; // Minimum point value post-convolution to accept as a seed (inclusive)
  int m_localMaxWindowSize = 0; // Only create roads from a local maximum, requires traceHits
  
  ///////////////////////////////////////////////////////////////////////
  // Convenience
  
  unsigned m_nLayers; 

  double m_step_x = 0; // step size of the bin boundaries in x
  double m_step_y = 0; // step size of the bin boundaries in y
  std::vector<double> m_bins_x; // size == m_imageSize_x + 1.
  std::vector<double> m_bins_y; // size == m_imageSize_y + 1
  // Bin boundaries, where m_bins_x[i] is the lower bound of bin i.
  // These are calculated from m_xMin/m_xMax
  
  ///////////////////////////////////////////////////////////////////////
  // Event Storage
  
  Image m_image;
  
  ///////////////////////////////////////////////////////////////////////
  // Core

  
  Image createLayerImage(unsigned layer, std::vector<const SimSpacePoint*> const & spacepoints) const;
  Image createImage(std::vector<const SimSpacePoint*> const & spacepoints) const;

  ///////////////////////////////////////////////////////////////////////
  // Helpers
  
  std::pair<unsigned, unsigned> yToXBins(size_t yBin_min, size_t yBin_max, double r, double phi, unsigned layer) const;
  unsigned getExtension(unsigned y, unsigned layer) const;
  bool passThreshold(Image const & image, unsigned x, unsigned y) const;
  void drawImage(Image const & image, std::string const & name);
};

}  // namespace ActsExamples




};



#endif // HTTHOUGHTRANSFORMTOOL_H
