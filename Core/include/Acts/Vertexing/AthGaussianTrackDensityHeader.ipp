/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRKVERTEXSEEDFINDERUTILIS_GAUSSIANTRACKDENSITY_H
#define TRKVERTEXSEEDFINDERUTILIS_GAUSSIANTRACKDENSITY_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "TrkVertexFitterInterfaces/IVertexTrackDensityEstimator.h"

#include <map>

namespace Trk
{

  class Track;
  class GaussianTrackDensity;

  /**
   @class GaussianTrackDensity

   Implementation of IVertexTrackDensityEstimator modeling reconstructed tracks as
   two-dimensional Gaussian distributions in (d0, z0) space and sampling
   the aggregate density distribution at user-requested points along the beam axis.

   @author Dave Casper <dcasper@uci.edu>

   ------------------------------------
   Changes:
  */

  class GaussianTrackDensity : public extends<AthAlgTool, IVertexTrackDensityEstimator>
  {
  public:
    /// Inherit constructor.
    using base_class::base_class;


    /**
     * @brief Find position of global maximum for density function.
     * @param vectorTrk List of input tracks.
     */
    virtual double
    globalMaximum (const std::vector<const Track*>& vectorTrk) const override;


    /**
     * @brief Find position of global maximum for density function.
     * @param vectorTrk List of input tracks.
     * @param density[out] Helper to hold density results.
     */
    virtual double
    globalMaximum (const std::vector<const Track*>& vectorTrk,
                   std::unique_ptr<ITrackDensity>& density) const override;


    /**
     * @brief Find position of global maximum for density function.
     * @param perigeeList List of input tracks.
     */
    virtual double
    globalMaximum (const std::vector<const TrackParameters*>& perigeeList) const override;


    /**
     * @brief Find position of global maximum for density function.
     * @param perigeeList List of input tracks.
     * @param density[out] Helper to hold density results.
     */
    virtual double
    globalMaximum (const std::vector<const TrackParameters*>& perigeeList,
                   std::unique_ptr<ITrackDensity>& density) const override;

    virtual
    std::pair<double,double> globalMaximumWithWidth (const std::vector<const TrackParameters*>& perigeeList/*,
                                                  std::unique_ptr<ITrackDensity>& density*/) const override;


  private:
    class TrackDensity;


    /**
     * @brief Find position of global maximum for density function.
     * @param pergigeeList List of input tracks.
     * @param density Helper density object.
     */
    double
    globalMaximumImpl (const std::vector<const TrackParameters*>& perigeeList,
                       TrackDensity& density) const;

    /**
     * @brief Find position of global maximum with Gaussian width for density function.
     * @param pergigeeList List of input tracks.
     * @param density Helper density object.
     */
    std::pair<double,double>
    globalMaximumWithWidthImpl (const std::vector<const TrackParameters*>& perigeeList,
                       TrackDensity& density) const;


    /**
     * @brief Add a set of tracks to a density object.
     * @param perigeeList Set of track parameters to add.
     * @param density Density object to which to add.
     */
    void addTracks(const std::vector<const TrackParameters*>& perigeeList,
                   TrackDensity& density) const;


    // functor to compare two Perigee values
    struct pred_perigee {
      bool operator()(const Perigee& left, const Perigee& right) const
      {
	return left.parameters()[Trk::z0] < right.parameters()[Trk::z0];
      }
    };

    struct TrackEntry
    {
      TrackEntry() { c_0 = 0; c_1 = 0; c_2 = 0; lowerBound = 0; upperBound = 0; }
      TrackEntry(double c0, double c1, double c2, double zMin, double zMax);
      TrackEntry(double zProbe);
      // Cached information for a single track
      double c_0;      // z-independent term in exponent
      double c_1;      // linear coefficient in exponent
      double c_2;      // quadratic coefficient in exponent
      double lowerBound;
      double upperBound;
    };

    // functor to compare two TrackEntry values based on their lower limits (low to high)
    struct pred_entry_by_min {
      bool operator()(const TrackEntry& left, const TrackEntry& right) const
      {
	return left.lowerBound < right.lowerBound;
      }
    };

    // functor to compare two TrackEntry values based on their upper limits (low to high)
    struct pred_entry_by_max {
      bool operator()(const TrackEntry& left, const TrackEntry& right) const
      {
	return left.upperBound < right.upperBound;
      }
    };

    typedef std::map< Perigee, 
                      GaussianTrackDensity::TrackEntry, 
                      GaussianTrackDensity::pred_perigee > trackMap;

    typedef std::map< Perigee, 
                      GaussianTrackDensity::TrackEntry, 
                      GaussianTrackDensity::pred_perigee >::const_iterator trackMapIterator;

    typedef std::map< GaussianTrackDensity::TrackEntry,
                      Perigee,
                      GaussianTrackDensity::pred_entry_by_max > lowerMap;

    typedef std::map< GaussianTrackDensity::TrackEntry,
                      Perigee,
                      GaussianTrackDensity::pred_entry_by_max >::const_iterator lowerMapIterator;

    typedef std::map< GaussianTrackDensity::TrackEntry,
                      Perigee,
                      GaussianTrackDensity::pred_entry_by_min > upperMap;

    typedef std::map< GaussianTrackDensity::TrackEntry,
                      Perigee,
                      GaussianTrackDensity::pred_entry_by_min >::const_iterator upperMapIterator;


    class TrackDensity : public ITrackDensity
    {
    public:
      TrackDensity (bool gaussStep) : m_gaussStep (gaussStep) {}
      virtual ~TrackDensity() = default;


      /**
       *  Evaluate the density function at the specified coordinate
       *  along the beamline.
       */
      virtual double trackDensity (double z) const override;


      /**
       *  Evaluate the density and its first two derivatives
       *  at the specified coordinate.
       */
      virtual void trackDensity (double z,
                                 double& density,
                                 double& firstDerivative,
                                 double& secondDerivative) const override;


      /**
       * @brief Return position of global maximum for density function.
       * @param msg Message stream.
       */
      double globalMaximum (MsgStream& msg) const;

      /**
       * @brief Return position of global maximum with Gaussian width for density function.
       * @param msg Message stream.
       */
      std::pair<double,double> globalMaximumWithWidth (MsgStream& msg) const;


      /**
       * @brief Add a track to the set being considered.
       * @param itrk Track parameters.
       * @param d0SignificanceCut Significance cut on d0.
       * @param z0SignificanceCut Significance cut on z0.
       */
      void addTrack (const Perigee& itrk,
                     const double d0SignificanceCut,
                     const double z0SignificanceCut);




    private:
      inline void updateMaximum(double trialZ, double trialValue, double secondDerivative, double& maxZ, double& maxValue, double& maxSecondDerivative) const
      { if (trialValue > maxValue) { maxZ = trialZ; maxValue = trialValue; maxSecondDerivative=secondDerivative;} }

      inline double stepSize(double y, double dy, double ddy) const
      { return ( m_gaussStep ? (y * dy)/(dy * dy - y * ddy) : -dy/ddy ); }

      bool m_gaussStep;

      //  Cache for track information
      trackMap m_trackMap;
      lowerMap m_lowerMap;
      upperMap m_upperMap;

      double m_maxRange = 0;
    };

    
    //  Cuts set by configurable properties
    
    //  Maximum allowed d0 significance to use (in sigma)
    Gaudi::Property<double> m_d0MaxSignificance { this, 
                                                  "MaxD0Significance", 
	                                          3.5, 
                                                  "Maximum radial impact parameter significance to use track" };

    // Tracks within this many sigma(z) are added to weight; increasing cut trades CPU efficiency for improved smoothness in tails
    Gaudi::Property<double> m_z0MaxSignificance { this,
	                                          "MaxZ0Significance",
	                                          12.0,
	                                          "Maximum longitudinal impact parameter significance to include track in weight"};

    // Assumed shape of density function near peak; Gaussian (true) or parabolic (false)
    Gaudi::Property<bool> m_gaussStep           { this,
	                                          "GaussianStep",
	                                          true,
	                                          "Peak search: True means assume Gaussian behavior, False means Newton/parabolic" };
                            
  };
}
#endif
