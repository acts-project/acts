///////////////////////////////////////////////////////////////////
// Layer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_LAYER_H
#define ACTS_DETECTOR_LAYER_H 1

// Core module
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometryObject.hpp"
#include "ACTS/Utilities/ApproachDescriptor.hpp"
#include "ACTS/Utilities/OverlapDescriptor.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Intersection.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"

namespace Acts {
  
  class Surface;
  class SurfaceMaterial;
  class MaterialProperties;
  class BinUtility;
  class Volume;
  class VolumeBounds;
  class TrackingVolume;
  class DetachedTrackingVolume;
  class ApproachDescriptor;
  class ICompatibilityEstimator;

  typedef ObjectIntersection<Surface> SurfaceIntersection;

  // master typedef
  class Layer;
  typedef std::shared_ptr<const Layer> LayerPtr;
  typedef std::pair<const Layer*, const Layer*> NextLayers;
  

  /**
     @enum LayerType

     For readability
  */
  enum LayerType { passive = 0,
                   active = 1 };

  /**
     @class Layer

     Base Class for a Detector Layer in the Tracking realm.
     An actual implemented Detector Layer inheriting from this base
     class has to inherit from a specific type of Surface as well.
     In addition, a Layer can carry:

     - a SurfaceArray of Surfaces holding the actual detector elements or subSurfaces.
     - SurfaceMaterial for Surface-based materialUpdates
     - an OverlapDescriptor (mainly used for blind extrapolation)
     - a pointer to the TrackingVolume (can only be set by such)
     - an active/passive code :
     0      - activ
     1      - passive
     [....] - other

     The search type for compatible surfaces on a layer is [ the higher the number, the faster ]:
     --------- Layer internal ------------------------------------------------               
     -1     - untested: provide all layer surfaces to the extrapolation engine 
                   - does not work with endSurface, will be increased to 0 if endSurface is given
                   - debug mode only !
      0     - test all on intersection and provide to the extrapolation engine
     --------- Overlap descriptor --------------------------------------------              
      1     - provide bin surface and registered neighbours and bin mates
                   - does not work with endSurface, will be increased to 2 if endSurface is given
      2     - as 1 but with intersection test @TODO compatibility test
      3     - provide bin surface and next bin surfaces (if differ)
                   - does not work with endSurface, will be increased to 4 if endSurface is given
      4     - as 3 but with intersection test @TODO compatibility test
      5     - whatever the overlap descriptor returns with this

     @author Andreas.Salzburger@cern.ch
  */

  class Layer : public virtual GeometryObject {

    /** Declare the TrackingVolume as a friend, to be able to register previous,
       next and set the enclosing TrackingVolume*/
    friend class TrackingVolume;

    /** Declare the DetachedTrackingVolume as a friend to be able to register it */
    friend class DetachedTrackingVolume;

  public:

    /** Clone at a with a shift - this is the only clone allowed */
    virtual LayerPtr cloneWithShift(const Transform3D& shift) const = 0;
    
    /** Destructor*/
    virtual ~Layer();

    /** Assignment operator - forbidden, layer assignment can be ambiguous */
    Layer& operator=(const Layer& lay) = delete;

    /** Return the entire SurfaceArray, returns a nullptr if no SurfaceArray */
    const SurfaceArray* surfaceArray() const;

    /** Transforms the layer into a Surface representation for extrapolation */
    virtual const Surface& surfaceRepresentation() const = 0;

    /** Return the Thickness of the Layer */
    double thickness() const;

    /** templated onLayer() method */
    template <class T> bool onLayer(const T& parameters, const BoundaryCheck& bcheck  = BoundaryCheck(true)) const;

    /** isOnLayer() method, using isOnSurface() with Layer specific tolerance */
    virtual bool isOnLayer(const Vector3D& gp, const BoundaryCheck& bcheck = BoundaryCheck(true)) const;

    /** getting the overlap descriptor */
    const OverlapDescriptor* overlapDescriptor() const;

    /** getting the approach descriptor */
    const ApproachDescriptor* approachDescriptor() const;

    /** Surface seen on approach - if not defined differently, it is the surfaceRepresentation() */
    virtual const SurfaceIntersection surfaceOnApproach(const Vector3D& pos,
                                                        const Vector3D& dir,
                                                        PropDirection pdir,
                                                        const BoundaryCheck& bcheck,
                                                        bool resolveSubSurfaces = false,
                                                        const ICompatibilityEstimator* ice = nullptr) const;

    /** get compatible surfaces starting from charged parameters
        returns back the compatible surfaces either with straight line estimation,
        or (@TODO later) with a compatiblityEstimator.
        - if start/end surface are given, surfaces are provided in between (start & end excluded)
        - the boolean indicates if the surfaces are direction ordered
     */
    virtual bool compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
			                        const TrackParameters& pars,
			                        PropDirection pdir,
			                        const BoundaryCheck& bcheck,
                                    bool collectSensitive, 
                                    bool collectPassive, 
                                    int searchType, 
			                        const Surface* startSurface = nullptr,
			                        const Surface* endSurface = nullptr,
			                        const ICompatibilityEstimator* ice = nullptr) const;

    /** get compatible surfaces starting from charged parameters
        returns back the compatible surfaces either with straight line estimation,
        or (@TODO later) with a compatiblityEstimator.                                                                                                                                                                                                                                      
        - if start/end surface are given, surfaces are provided in between (start & end excluded)
        - the boolean indicates if the surfaces are direction ordered
     */
    virtual bool compatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
    			                    const NeutralParameters& pars,
    			                    PropDirection pdir,
    			                    const BoundaryCheck& bcheck,
                                    bool collectSensitive, 
                                    bool collectPassive, 
                                    int searchType, 
    			                    const Surface* startSurface = nullptr,
    			                    const Surface* endSurface = nullptr,
    			                    const ICompatibilityEstimator* ice = nullptr) const;
                                      

    /** Has sub-structure method:
        - sub-structure depending on :
          (a) only when required to resolve sub surfaces for sensitive hits
          (b) also material is ordered with sub structure */
    virtual bool hasSubStructure(bool resolveSensitive=false) const;

    /** Boolean check method if layer has material:
       - checks if any of the layer surfaces has material:
       - can be approach surfaces or layer surface */
    virtual bool hasMaterial() const;

    /** Boolean check method if layer has sensitive surfaces */
    virtual bool hasSensitive() const;

    const Layer* nextLayer(const Vector3D& pos, const Vector3D& dir) const;

    /** get the confining TrackingVolume */
    const TrackingVolume* enclosingTrackingVolume() const;

    /** get the confining DetachedTrackingVolume */
    const DetachedTrackingVolume* enclosingDetachedTrackingVolume() const;

    /** register Volume associated to the layer - if you want to do that by hand*/
    void registerRepresentingVolume(const AbstractVolume *theVol) const;

    /** get the Volume associated to the layer */
    const AbstractVolume* representingVolume() const;

  protected:
    /** Default Constructor*/
    Layer();

    /** Constructor with pointer to SurfaceArray (passing ownership) */
    Layer(SurfaceArray* surfaceArray,
          double thickness = 0.,
          OverlapDescriptor* od = nullptr,
          ApproachDescriptor* ad = nullptr,
          int ltype=int(passive));

    /** get compatible surfaces starting from charged parameters - forward call from explicit methods */
    template <class T>  bool getCompatibleSurfaces(std::vector<SurfaceIntersection>& cSurfaces,
                                                   const T& pars,
                                                   PropDirection pdir,
                                                   const BoundaryCheck& bcheck,
                                                   bool collectSensitive, 
                                                   bool collectPassive, 
                                                   int searchType, 
                                                   const Surface* startSurface = nullptr,
                                                   const Surface* endSurface = nullptr,
                                                   const ICompatibilityEstimator* ice = nullptr) const;
                                                     
    /** test compatible surface - checking directly for intersection & collection */
    void testCompatibleSurface(std::vector<SurfaceIntersection>& cSurfaces,
                               const Surface& surface,
   	                           const Vector3D& pos,
                               const Vector3D& dir,
   	                           PropDirection pdir,
   	                           const BoundaryCheck& bcheck,
                               double maxPathLength, 
                               bool collectSensitive, 
                               bool collectPassive, 
                               bool intersectionTest,
   	                           const Surface* startSurface = nullptr,
   	                           const Surface* endSurface = nullptr,
   	                           const ICompatibilityEstimator* ice = nullptr) const;
   
                                                     

    /** private method to set enclosing TrackingVolume, called by friend class only
        optionally, the layer can be resized to the dimensions of the TrackingVolume
        - Bounds of the Surface are resized
        - MaterialProperties dimensions are resized
        - SubSurface array boundaries are NOT resized
    */
    void encloseTrackingVolume(const TrackingVolume& tvol) const;

    /** private method to set the enclosed detached TV,
       called by friend class only */
    void encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol) const;


    mutable NextLayers                              m_nextLayers;                       //!< the previous Layer according to BinGenUtils
    mutable const BinUtility*                       m_nextLayerUtility;                 //!< the bin utility to find the next layer

    SurfaceArray*                                   m_surfaceArray;                     //!< SurfaceArray on this layer Surface
    double                                          m_layerThickness;                   //!< thickness of the Layer

    OverlapDescriptor*                              m_overlapDescriptor;                //!< descriptor for overlap/next surface
    mutable ApproachDescriptor*                     m_approachDescriptor;               //!< surface for approaching
    mutable const TrackingVolume*                   m_enclosingTrackingVolume;          //!< Enclosing TrackingVolume
    mutable const DetachedTrackingVolume*           m_enclosingDetachedTrackingVolume;  //!< Enclosing DetachedTrackingVolume

    mutable const AbstractVolume*                   m_representingVolume;               //!< Representing Volume - the surfaces of that can be used as

    int                                             m_layerType;                        //!< make a passive/active division

  };

  inline const SurfaceArray* Layer::surfaceArray() const
    { return m_surfaceArray; }

  inline double Layer::thickness() const
    { return m_layerThickness; }
  
  inline const OverlapDescriptor* Layer::overlapDescriptor() const
    { return m_overlapDescriptor; }

  inline const TrackingVolume* Layer::enclosingTrackingVolume() const
    { return m_enclosingTrackingVolume; }

  inline void Layer::encloseTrackingVolume(const TrackingVolume& tvol) const
    { m_enclosingTrackingVolume = &(tvol); }

  inline const DetachedTrackingVolume* Layer::enclosingDetachedTrackingVolume() const
    { return m_enclosingDetachedTrackingVolume; }

  inline void Layer::encloseDetachedTrackingVolume(const DetachedTrackingVolume& tvol) const
    { m_enclosingDetachedTrackingVolume = &(tvol); }

  inline const AbstractVolume* Layer::representingVolume() const
    { return m_representingVolume; }

  inline const Layer* Layer::nextLayer(const Vector3D& gp, const Vector3D& mom) const {
     // no binutility -> no chance to find out the direction
     if (!m_nextLayerUtility) return nullptr;
     return (m_nextLayerUtility->nextDirection(gp, mom) < 0) ? m_nextLayers.first : m_nextLayers.second;
  }

  inline void Layer::registerRepresentingVolume(const AbstractVolume *theVol) const
    { delete m_representingVolume; m_representingVolume = theVol; }

  #include "ACTS/Layers/detail/Layer.icc"

  /** Layers are constructedd with shared_ptr factories, hence the layer array is describes as: */
  typedef BinnedArray< LayerPtr > LayerArray;
  
} // end of namespace

#endif // ACTS_DETECTOR_LAYER_H

