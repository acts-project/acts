/*
  Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TRKVERTEXSEEDFINDERTOOLS_TrackDensitySeedFinder_H
#define TRKVERTEXSEEDFINDERTOOLS_TrackDensitySeedFinder_H

#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "TrkVertexFitterInterfaces/IVertexAnalyticSeedFinder.h"
#include "TrkVertexFitterInterfaces/IVertexTrackDensityEstimator.h"

namespace Trk
{

  class Track;
  class IMode1dFinder;
  class ITrackToVertexIPEstimator;

  // @author D. Casper
  //
  // @ATLAS software
  //
  // This class implements a seed finder for the primary vertex finding 
  // which is based on the use of a density function along the beam line 
  // 
  // -------------------------------------------

  class TrackDensitySeedFinder : public extends<AthAlgTool, IVertexAnalyticSeedFinder>
  {
  public:
    // Standard Gaudi constructor.
    TrackDensitySeedFinder (const std::string& t,
                            const std::string& n,
                            const IInterface*  p);


    virtual ~TrackDensitySeedFinder();


    virtual StatusCode initialize() override;
    virtual StatusCode finalize() override;


    using IVertexSeedFinder::findSeed;

    /**
     *  Finds a linearization point out of a vector of tracks
     *  and returns it as an Amg::Vector3D object. If you want an 
     *  additional constraint can be taken into account.
     */
    virtual Amg::Vector3D
    findSeed (const std::vector<const Trk::Track*> & vectorTrk,
              const xAOD::Vertex * constraint=0) const override;
    

    /** 
     * Finds a linearization point out of a vector of TrackParameters
     *  and returns it as an Amg::Vector3D object. If you want an 
     * additional constraint can be taken into account.
     */
    virtual Amg::Vector3D
    findSeed(const std::vector<const Trk::TrackParameters*> & perigeeList,
             const xAOD::Vertex * constraint=0) const override;

    virtual std::pair<Amg::Vector3D,Amg::MatrixX>  findAnalyticSeed (const std::vector<const Trk::TrackParameters*>& perigeeList,
    														   const xAOD::Vertex * constraint=0) const override;


    /**
     * Finds full vector of linearization points from a vector of tracks
     *  and returns it as an Amg::Vector3D object.  Intended for seed finders that produce all at once.
     *  If you want an additional constraint can be taken into account.
     */
    virtual std::vector<Amg::Vector3D>
    findMultiSeeds (const std::vector<const Trk::Track*>& vectorTrk,
                    const xAOD::Vertex * constraint=0) const override;


    /**
     * Finds full vector of linearization points from a vector
     * of TrackParameters and returns it as an Amg::Vector3D object.
     * Intended for seed finders that produce all at once.
     * If you want an additional constraint can be taken into account.
     */
    virtual std::vector<Amg::Vector3D>
    findMultiSeeds (const std::vector<const Trk::TrackParameters*>& perigeeList,
                    const xAOD::Vertex * constraint=0) const override;


  private:
    ToolHandle< IVertexTrackDensityEstimator > m_densityEstimator { this, 
	                                                            "DensityEstimator", 
	                                                            "Trk::GaussianTrackDensity",
                                                                    "Track density tool" } ;

  };
}
#endif
