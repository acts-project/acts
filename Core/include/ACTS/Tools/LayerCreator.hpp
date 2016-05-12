///////////////////////////////////////////////////////////////////
// LayerCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_LAYERCREATOR_H
#define ACTS_GEOMETRYTOOLS_LAYERCREATOR_H 1

// Geometry module
#include "ACTS/Tools/ILayerCreator.hpp"
#include "ACTS/Tools/ISurfaceArrayCreator.hpp"

#ifndef ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define ACTS_LAYERCREATOR_TAKESMALLERBIGGER
#define takeSmaller(current,test) current = current < test ? current : test
#define takeBigger(current,test)  current = current > test ? current : test
#define takeSmallerBigger(cSmallest, cBiggest, test) takeSmaller(cSmallest, test); takeBigger(cBiggest, test)
#endif

namespace Acts {
    
    class ISurfaceArrayCreator;
    
    /** @class LayerCreator
     
     The LayerCreator is able to build cylinde, disc layers or plane layers from detector elements
     
     @author Andreas.Salzburger@cern.ch
     */
    
    class LayerCreator : public ILayerCreator {
        
    public:
        /** @struct Config 
            Configuration for the LayerCreator */
        struct Config {
        
            std::shared_ptr<ISurfaceArrayCreator>           surfaceArrayCreator; //!< geometry tool binning the surfaces into arrays
            
            Config() :
                surfaceArrayCreator(nullptr)
            {}
        };
        /** constructor */
        LayerCreator(const Config& lcConfig);
        
        /** destructor */
        ~LayerCreator() = default;

        /** ILayerCreator interface method - returning a cylindrical layer */
        LayerPtr cylinderLayer(const std::vector<const Surface*>& surfaces,
                               double envelopeR, double evelopeZ,
                               size_t binsPhi, size_t binsZ) const override; 
      
        /** ILayerCreator interface method - returning a disc layer */
        LayerPtr discLayer(const std::vector<const Surface*>& surfaces,
                           double envelopeMinR, double envelopeMaxR, double envelopeZ,
                           size_t binsR, size_t binsZ,
                           const std::vector<double>& rBoundaries = {}) const override; 
      
        /** ILayerCreator interface method - returning a plane layer */
        LayerPtr planeLayer(const std::vector<const Surface*>& surfaces,
                            double envelopeXY, double envelopeZ,
                            size_t binMultiplierX, size_t binMultiplierY) const override;
        
        /* set the configuration object**/
        void setConfiguration(const Config& lcConfig);
        /** access th configuration object */
        Config getConfiguration() const;
    
    private:
        
        /** method to get the global extends in space for the module */
        void moduleExtend(const Surface& sf,
                          double& minR,   double& maxR, 
                          double& minPhi, double& maxPhi, 
                          double& minZ,   double& maxZ) const;        
        
        /** calculates the closest radial distance of a line */
        double radialDistance(const Vector3D& pos1, const Vector3D& pos2) const;
        
        /** configuration object */
        Config                                      m_config;

    };
    
    inline LayerCreator::Config LayerCreator::getConfiguration() const
        { return m_config; }
    
} //end of namespace

#endif // ACTS_GEOMETRYTOOLS_LAYERCREATOR_H
