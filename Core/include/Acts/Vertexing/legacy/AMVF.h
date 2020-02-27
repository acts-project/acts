/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @class InDet::InDetAdaptiveMultiPriVxFinderTool
 *
 * @author Giacinto Piacquadio (Freiburg University)
 * 
 * (giacinto.piacquadio@physik.uni-freiburg.de)
 * 
 * This class provides an implementation for a primary 
 * vertex finding tool, which uses the \<i\>Adaptive MultiVertex\</i\> 
 * Fitter to solve the problem of finding the multiple 
 * interaction vertices and to find out the best possible 
 * assignment of the track to the vertices.
 *
 * The steps done are the following:
 * - the selection cuts are applied
 * 
 * then iteratively:
 * - a new vertex is seeded with the remaining tracks
 *    (the seed finder is used)
 * - all the tracks whose Z at PCA is closer to the seeded 
 *    vertex by more than TracksMaxZinterval (by JobOption), 
 *    are added to the new vertex candidate
 * -  the new vertex candidate is added on top of the previous fit and 
 *    the AdaptiveMultiVertexFitter is used to fit all them 
 *    together (more information in the \<br\>\<i\>TrkVertexFitters\</br\>\</i\>
 *    package).
 * -  the tracks already used are removed from the tracks 
 *    from which the next seed would be obtained and if there 
 *    are more than 2 left, a new iteration is started.
 *
 * when no more than 2 seeding tracks are left:
 * - a vector of MVFVxCandidate is provided as result and 
 *   according to the selection type, the order in which it is 
 *   provided represents how high the probability of that 
 *   particular vertex to come from the primary vertex is. 
 *
 * In general the first VxCandidate* in the collection is 
 * the one with highest sqrt(N_tracks)*Sum Pt_track^2. This 
 * is the case if the selectiontype in the jobOptions is set 
 * to 0 (default).
 *
 *
 * This finder is particularly suited for the high luminosities 
 * scenarios which will came up at LHC.
 *
 * ------------------------------------------------------------
 * Changes:
 *
 * David Shope <david.richard.shope@cern.ch> (2016-04-19)
 *
 * EDM Migration to xAOD - from Trk::VxCandidate to xAOD::Vertex, 
 *                         from Trk::RecVertex   to xAOD::Vertex,
 *                         from Trk::Vertex      to Amg::Vector3D
 *
 * Also, VxMultiVertex EDM has been migrated to the following:
 *
 *   Trk::MvfFitInfo  
 *     constraintVertex     now uses xAOD::Vertex
 *     seedVertex           now uses Amg::Vector3D
 *     linearizationVertex  now uses Amg::Vector3D
 *
 *   Trk::TrackToVtxLink
 *     Vertex objects stored using this class are now xAOD::Vertex
 *
 * Instead of using the MVFVxCandidate class, xAOD::Vertex is employed by decorating it
 * with the multi-vertex information:
 *
 *   bool                              isInitialized
 *   MvfFitInfo*                       MvfFitInfo
 *   std::Vector\<VxTrackAtVertex*\>     VTAV
 *
 *   This last decoration is needed in order to be able to use MVFVxTrackAtVertex objects
 *   which have the additional information of xAOD::Vertices associated to the track
 *   and (through polymorphism) to still be able to pass them to the KalmanVertexUpdator as
 *   VxTracksAtVertex objects.
 */

#ifndef INDETPRIVXFINDERTOOL_INDETADAPTIVEMULTIPRIVXFINDERTOOL_H
#define INDETPRIVXFINDERTOOL_INDETADAPTIVEMULTIPRIVXFINDERTOOL_H

#include "InDetRecToolInterfaces/IVertexFinder.h"
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"
#include "TrkTrack/TrackCollection.h" // type def ...
#include "TrkParticleBase/TrackParticleBaseCollection.h"
#include "TrkParameters/TrackParameters.h"

/**
 * Forward declarations 
 */
#include "xAODTracking/VertexFwd.h"
#include "xAODTracking/TrackParticleFwd.h"
#include "xAODTracking/VertexContainerFwd.h"
#include "xAODTracking/TrackParticleContainerFwd.h"
#include "BeamSpotConditionsData/BeamSpotData.h"
#include "TrkVertexFitterInterfaces/ITrackToVertexIPEstimator.h"
class TrackToVtxLinkContainer;
class NN;


namespace Trk
{
  class IVertexAnalyticSeedFinder;
  class AdaptiveMultiVertexFitter;
  class Track;
  class ITrackLink;
  class TrkQuality;
  class IVxCandidateXAODVertex;
}


namespace InDet
{
  class IInDetTrackSelectionTool;

 class InDetAdaptiveMultiPriVxFinderTool : public AthAlgTool, virtual public IVertexFinder
 {
 
  public:
    
    InDetAdaptiveMultiPriVxFinderTool(const std::string& t, const std::string& n, const IInterface*  p);
    virtual ~InDetAdaptiveMultiPriVxFinderTool();
   
    StatusCode initialize();

    /**
     * The MultiVertexFinding is performed.
     *
     * Input is the Track Collection. Output is the VertexContainer 
     * with a list of fitted vertices, according to the probability 
     * of being the primary interaction point.
     * 
     * Description of the finder is provided as doxygen info on the class.
     *
     */

    std::pair<xAOD::VertexContainer*, xAOD::VertexAuxContainer*> findVertex(const TrackCollection* trackTES);
    std::pair<xAOD::VertexContainer*, xAOD::VertexAuxContainer*> findVertex(const Trk::TrackParticleBaseCollection* trackTES);
    std::pair<xAOD::VertexContainer*, xAOD::VertexAuxContainer*> findVertex(const xAOD::TrackParticleContainer* trackParticles);

    StatusCode finalize();
    
  private:

    std::pair<xAOD::VertexContainer*, xAOD::VertexAuxContainer*> findVertex(const std::vector<const Trk::ITrackLink*> & trackVector);
  
    void SGError(std::string errService);
    virtual void printParameterSettings();

    ToolHandle< Trk::AdaptiveMultiVertexFitter > m_MultiVertexFitter;
    ToolHandle< Trk::IVertexAnalyticSeedFinder > m_analyticSeedFinder;
    ToolHandle< InDet::IInDetTrackSelectionTool > m_trkFilter;

    SG::ReadCondHandleKey<InDet::BeamSpotData> m_beamSpotKey { this, "BeamSpotKey", "BeamSpotData", "SG key for beam spot" };

    /** Define a beam constraint for the fit */
    bool m_useBeamConstraint; //!<  Use a vertex/beam constraint
    

    /**
     * When adding a new vertex to the multi vertex fit,
     * only the tracks whose Z at PCA is closer 
     * to the seeded by more than this TracksMaxZinterval 
     * value are added to this new vertex.
     *
     * Default is 4 mm. If you cut too hard, you cut out 
     * the good cases where the seed finder is not 
     * reliable, but the fit would be still able to converge 
     * towards the right vertex. If you cut too soft, you 
     * consider a lot of tracks which just slow down the fit.
     */

    double m_TracksMaxZinterval;

    /**
     * After having added one vertex to the fit and having 
     * performed the MultiVertex fit again, all the tracks 
     * which are compatible to the new vertex by more than 
     * this maxVertexChi2 (in units of chi2) value are eliminated from the 
     * tracks from which still to seed the next vertex.
     *
     */

    double m_maxVertexChi2;

    /**
     * As a default the realMultiVertex should stay to false (because this is very well tested).
     *
     * If switched to true, all the tracks are considered to be added to the new vertex 
     * after this new one is seeded, and not only the ones which are considered as outliers 
     * of previous fitted vertices.
     *
     * The presence of a core of tracks the previous vertices are as attached to stabilizes 
     * the fit quite drastically. In case of luminosities higher than the low lumi scenario,
     * one should probably to try to switch this on, or, if this doesn't work, decrease the 
     * maxVertexChi2 and the cleaningZinterval to lower values.
     */
    
    bool m_realMultiVertex;


    /*
     * Decides if you want to use the vtxCompatibility() of the track (set to true) or 
     * the chi2() (set to false) as an estimate for a track being an outlier or not.
     * The vtxCompatibility() is the default. In case the track refitting 
     * is switched on in the AdaptiveMultiVertex fitter, you may want to 
     * use the refutted chi2().
     *
     */

    bool m_useFastCompatibility;

    /*
     * Selection of the most probable primary interaction vertex is done accordingly to:
     * - selectiontype is 0:  just sqrt(N_tracks)*Sum_track p_t_track^2
     * - selectiontype is 1:  Neural Network is used, trained on WH(120)
     */

    int m_selectiontype;

    /*
     * During the estimation of probability of vertex candidate to be the primary interaction 
     * vertex, only all the tracks which have chi2 in the vertex fit higher than this value 
     * are used for the sum of p_t^2 or as input for the Neural Network.
     */

    double m_finalCutMaxVertexChi2;

 
    /*
     * Maximum significance on the distance between two vertices 
     * to allow merging of two vertices.
     *
     */

    double m_cutVertexDependence;
    

    /*
     * Has to be setup equal to the minimum weight set in the fitter.
     *
     * In the fitting, when a track has a weight lower than this value,
     * the track is not updated during that iteration.
     */

    double m_minweight;


    /*
    * Impact parameter estimator used to calculate significance
    */
    ToolHandle< Trk::ITrackToVertexIPEstimator > m_ipEstimator { "Trk::TrackToVertexIPEstimator" };
    
    /*
     * Maximum amount of iterations allowed for vertex finding.
     * 
     * The more vertices you have in the event, the more iterations you have to 
     * allow (safe factor: number of expected vertices * 10)
     *
     */

    double m_maxIterations;

    NN* m_testingclass;

   /*
    * Fit also single track vertices
    * (could be usefull for example for H-> gamma gamma)\
    *
    */

   bool m_addSingleTrackVertices;

   bool m_do3dSplitting;
   
   double m_zBfieldApprox;

   double m_maximumVertexContamination;

    /*
    * Maximum allowed significance of track position to vertex seed
    */
    double m_tracksMaxSignificance ;

    /*
    * Toggle vertex seed constraint on/off
    */
    bool m_useSeedConstraint ;


   struct CompareTheTwoVertices {
     bool operator()( xAOD::Vertex* const & first, xAOD::Vertex* const & second);
   };

   /**												     
    * Internal method to estimate the probability to be signal vertex of a certain vertex candidate. 
    */												     

   double estimateSignalCompatibility(xAOD::Vertex *myxAODVertex);				     

   /**
    * Estimate DeltaZ given a certain track parameters and beam spot center position
    * ONLY TEMPORARY 15-08-2009: common tool needed to collect this method
    */
   
   double estimateDeltaZ(const Trk::TrackParameters& myPerigee, const Amg::Vector3D& myTransvVertex);

   /** 
   * copying from the guassian density alg
   */
   double ipSignificance(const Trk::TrackParameters* params, const Amg::Vector3D * vertex) const;

   /**
    * Clean decorator data from a vertex candidate (to avoid memory leaks) and then delete it and set to zero
    */
   
   void releaseCandidate(xAOD::Vertex*& candidate);
   

 };//end of class
}//end of namespace definitions
#endif
