///////////////////////////////////////////////////////////////////
// SurfaceArrayCreator.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Tools/SurfaceArrayCreator.hpp"
#include "ACTS/Utilities/BinStepping.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/BinnedArray2D.hpp"
#include "ACTS/Utilities/BinnedArray1D.hpp"
#include "ACTS/Surfaces/Surface.hpp"
// Core module
#include "ACTS/Utilities/Definitions.hpp"

// constructor
Acts::SurfaceArrayCreator::SurfaceArrayCreator()
{}

/** SurfaceArrayCreator interface method - create an array in a cylinder, binned in phi, z */
std::unique_ptr<Acts::SurfaceArray> Acts::SurfaceArrayCreator::surfaceArrayOnCylinder(const std::vector<const Acts::Surface*>& surfaces,
                                                                      double R, double minPhi, double maxPhi, double halfZ, 
                                                                      size_t binsPhi, size_t binsZ, 
                                                                      std::shared_ptr<Acts::Transform3D> transform) const
{

    //MSG_DEBUG("Creating a SurfaceArray on a cylinder with bins in phi x z = " << binsPhi << " x " << binsZ );

    // create the (plain) binUtility - with the transform
    BinUtility* arrayUtility = new BinUtility(binsPhi, minPhi, maxPhi, closed, binPhi, transform);
    (*arrayUtility) += BinUtility(binsZ, -halfZ, halfZ, open, binZ);

    // the z step and the phi step
    double zStep   = (2*halfZ/binsZ);
    double phiStep = (maxPhi-minPhi)/binsPhi;
    
    // prepare the surface system in phi x z  
    std::vector< std::vector< std::pair< SurfacePosition, Vector3D > > > phizSystem;
    phizSystem.reserve(binsZ);
    for (size_t iZ = 0; iZ < binsZ; ++iZ){
        // the current z value
        double currentZ = -halfZ + (iZ+0.5)*zStep;
        // the current phi row 
        std::vector< std::pair< SurfacePosition, Vector3D > > phiSystem;
        phiSystem.reserve(binsPhi);
        for (size_t iPhi = 0; iPhi < binsPhi; ++iPhi){
            // the current phi value
            double currentPhi = minPhi + (iPhi+0.5)*phiStep;
            // the bin position 
            Vector3D binPosition = transform ? ((*transform)*Vector3D(R*cos(currentPhi),R*sin(currentPhi),currentZ)) :
                                               Vector3D(R*cos(currentPhi),R*sin(currentPhi),currentZ);
            // the bin direction 
            Vector3D binDirection = transform ? ((transform->linear())*Vector3D(cos(currentPhi),sin(currentPhi),0.)) :
                                               Vector3D(cos(currentPhi),sin(currentPhi),0.);
            // push it in
            phiSystem.push_back( std::pair< SurfacePosition, Vector3D >(SurfacePosition(nullptr, binPosition), binDirection) );                                                                      
            
        }
        phizSystem.push_back(phiSystem);
    }
    
    // create and complete 
    std::vector<SurfacePosition> sVector;
    // complete
    completeBinning(surfaces, *arrayUtility, sVector, phizSystem);
    // create the surfaceArray
    auto sArray = std::make_unique<BinnedArray2D<const Surface*>>(sVector,arrayUtility);
    // register the neighbours  
    registerNeighboursGrid(sArray->arrayObjectsOrdered(), false, true);
    // return the surface array
    return std::move(sArray);
} 

/** SurfaceArrayCreator interface method - create an array on a disc, binned in r, phi */
std::unique_ptr<Acts::SurfaceArray> Acts::SurfaceArrayCreator::surfaceArrayOnDisc(const std::vector<const Acts::Surface*>& surfaces,
                                                                  double minR, double maxR, double minPhi, double maxPhi,
                                                                  size_t binsR, size_t binsPhi,
                                                                  const std::vector<double>& rBoundaries,
                                                                  std::shared_ptr<Acts::Transform3D> transform) const
{

    //MSG_DEBUG("Creating a SurfaceArray on a disc with bins in r x phi = " << binsR << " x " << binsPhi );
    
    // the z step and the phi step
    double phiStep = (maxPhi-minPhi)/(binsPhi-1);
    std::vector<SurfacePosition> sVector;
    
    // create the (plain) binUtility - with the transform
    BinUtility* arrayUtility = nullptr;
    if (binsR == 1){
        //MSG_DEBUG("Only one single ring is R is present - 1D surface array to be created." );
        // reserve the right amoung
        sVector.reserve(binsPhi);
        // solve this at once 
        arrayUtility  = new BinUtility(binsPhi, minPhi, maxPhi, closed, binPhi, transform);
        // fill the surfaces
        for (auto& surface : surfaces){
            // fill it into the surface,position for further registration
            sVector.push_back(SurfacePosition(surface,surface->binningPosition(binPhi)));
        }
        // create the surface array
        auto sArray = std::make_unique<BinnedArray2D<const Surface*>>(sVector,arrayUtility);
        // register the neighbours
        // copy to get const away, 
        // - but take them from the surface array (and not the input vector) because like this they are bin ordered
        std::vector<const Surface*> arraySurfaces;
        arraySurfaces.insert(arraySurfaces.begin(),sArray->arrayObjects().begin(),sArray->arrayObjects().end());
        std::vector< std::vector< const Surface*> > arraySystem = { arraySurfaces };
        // prepared to run the neighbour registration now
        registerNeighboursGrid(arraySystem, false, true);
        // now return
        return std::move(sArray);
    } 
    // more complicated binning 2D 
    double rStep   = ((maxR-minR)/(binsR));
    // 2D binning
    arrayUtility = new BinUtility(binsR, minR, maxR, open, binR, transform);
    (*arrayUtility) += BinUtility(binsPhi, minPhi, maxPhi, closed, binPhi);
    
    // prepare the surface system in r x phi  
    std::vector< std::vector< std::pair< SurfacePosition, Vector3D > > > rphiSystem;
    // the bin direction / @TODO should actually be +1, -1, gets important for intersection sorting method
    Vector3D binDirection(0.,0.,1);    
    rphiSystem.reserve(binsPhi);
    // loop of iR and iPhi bins and fill the order positions
    for (size_t iPhi = 0; iPhi < binsPhi; ++iPhi){
        // the current phi value
        double currentPhi = minPhi + (iPhi+0.5)*phiStep;
        // the current phi row 
        std::vector< std::pair< SurfacePosition, Vector3D > > rSystem;
        rSystem.reserve(binsR);
        for (size_t iR = 0; iR < binsR; ++iR){
            // the current R value
            double currentR = minR + (iR+0.5)*rStep;
            // the bin position 
            Vector3D binPosition = transform ? ((*transform)*Vector3D(currentR*cos(currentPhi),currentR*sin(currentPhi),0.)) :
                                               Vector3D(currentR*cos(currentPhi),currentR*sin(currentPhi),0.);
            // push it in
            rSystem.push_back( std::pair< SurfacePosition, Vector3D >(SurfacePosition(nullptr, binPosition), binDirection) );                                                                      
        
        }
        rphiSystem.push_back(rSystem);
    }    
    // create and complete 
    completeBinning(surfaces, *arrayUtility, sVector, rphiSystem);
    // create the surfaceArray
    auto sArray = std::make_unique<BinnedArray2D<const Surface*>>(sVector,arrayUtility);
    // register the neighbours  
    registerNeighboursGrid(sArray->arrayObjectsOrdered(), false, true);
    // return the surface array
    return std::move(sArray);
}

/** SurfaceArrayCreator interface method - create an array on a plane */
std::unique_ptr<Acts::SurfaceArray> Acts::SurfaceArrayCreator::surfaceArrayOnPlane(const std::vector<const Acts::Surface*>& /*surfaces*/,
                                                                   double /*halflengthX*/, double /*halflengthY*/, 
                                                                   size_t /*binsX*/, size_t /*binsY*/,
                                                                   std::shared_ptr<Acts::Transform3D> /*transform*/) const
{
    //!< @TODO implement - take from ATLAS complex TRT builder
    return nullptr;
}


void Acts::SurfaceArrayCreator::completeBinning(const std::vector<const Surface*>& surfaces, 
                                                const BinUtility& binUtility,
                                                std::vector<SurfacePosition>& sVector, 
                                                std::vector< std::vector< SurfacePositionDirection > >& binSystem) const 
{
    
    //MSG_DEBUG("Complete binning by filling closest neighbour surfaces in potentially empty bins." );
    
    // get the number of bins 
    size_t bins0 = binSystem.at(0).size();
    size_t bins1 = binSystem.size();

    //MSG_VERBOSE("Prefilling a bin system with [ " << bins0 << " x " << bins1 << " ].");

    // prefill the easy ones
    for (auto& surface : surfaces){
        // skip if not there
        if (!surface) continue;
        // calculate the bin
        size_t bin0 = binUtility.bin(surface->center(),0);
        size_t bin1 = binUtility.bin(surface->center(),1);
        //MSG_VERBOSE("- estimated bin [ " << bin0 << " x " << bin1 << " ] for surface " << surface);
        // fill it
        binSystem.at(bin1).at(bin0).first.first = surface;
    }
    
    size_t completedBins = 0;
        
    // now complete the system - register the neighbors if necessary
    for (size_t ibin1=0; ibin1 < bins1; ++ibin1){
        for (size_t ibin0 = 0; ibin0 < bins0; ++ibin0){
            // get the current surface 
            const Surface* binSurface = binSystem.at(ibin1).at(ibin0).first.first;
            // we are done when we have a surface registerd
            if (binSurface) continue;
            // binPosition
            Vector3D binPosition = binSystem.at(ibin1).at(ibin0).first.second;
            double surfaceDist = 10e10;   
            // counter 
            ++completedBins;
            // brute force method 
            // - find the closest surface 
            // @TODO try to add intersection test
            for (auto& surface : surfaces){
                // skip if not here
                if (!surface) continue;
                // recalculate distance
                double testDist = (binPosition-surface->center()).mag();
                if (testDist < surfaceDist){
                    binSystem.at(ibin1).at(ibin0).first.first = surface;
                    surfaceDist = testDist;
                }
            }
        }
    }
    //MSG_VERBOSE("Number of empty bins that were filled with neighbours: " << completedBins );
    
    // stream out the system for Binned Array creation
    sVector.reserve(bins0*bins1);
    for (auto& outer : binSystem)
        for (auto& inner : outer){
            //size_t bin0 = binUtility.bin(inner.first.second,0);
            //size_t bin1 = binUtility.bin(inner.first.second,1);
            //MSG_VERBOSE("- bin [ " << bin0 << " x " << bin1 << " ] holds surface " << inner.first.first);
            sVector.push_back(SurfacePosition(inner.first.first,inner.first.second));
       }
}

void Acts::SurfaceArrayCreator::registerNeighboursGrid(const std::vector< std::vector < const Surface* > >& surfaceArrayObjects, bool open0, bool open1) const {
    
    // get the number of bins
    size_t bins0 = surfaceArrayObjects.at(0).size();
    size_t bins1 = surfaceArrayObjects.size();
    
    size_t emptybins = 0;
    
    //MSG_DEBUG("Registering the neigbours for " << bins0 << " x " << bins1 << " bins.");
    
    // neighbour registration
    for (size_t i1 = 0; i1 < bins1; ++i1){
        for (size_t i0 = 0; i0 < bins0; ++i0){
            // the current one
            const Surface* surface  = surfaceArrayObjects.at(i1).at(i0);
            // leavit if there's no chance to do anything
            if (!surface || !surface->associatedDetectorElement()){
                ++emptybins;
                //MSG_WARNING(" - empty bin detected at [" << i0 << "][" << i1 << "]");
                continue;
            }    
            //MSG_VERBOSE("Processing neighbour registration for bin [" << i0 << "][" << i1 << "] - surface is " << surface);
            
            // the surface is defined
            const Surface* nsurface = surface;
            size_t p0 = i0;
            size_t n0 = i0;
            size_t p1 = i1;
            size_t n1 = i1;
            // find the previous neighbours in Loc0
            while (decrement(p0,bins0,open0)){
                nsurface = surfaceArrayObjects.at(i1).at(p0);
                //MSG_VERBOSE("                       - decrement to bin [" << p0 << "][" << i1 << "] - surface is " << nsurface);
                if (nsurface && nsurface != surface) break;
            }       
            // find the next neighbour in Loc0
            nsurface = surface;
            while (increment(n0,bins0,open0) && surface == nsurface){
                nsurface = surfaceArrayObjects.at(i1).at(n0);
                //MSG_VERBOSE("                       - increment to bin [" << n0 << "][" << i1 << "] - surface is " << nsurface);
                if (nsurface && nsurface != surface) break;
            }
            // find the previous neighbour in Loc1
            nsurface = surface;
            while (decrement(p1,bins1,open1) && surface == nsurface){
                nsurface = surfaceArrayObjects.at(p1).at(i0);
                //MSG_VERBOSE("                       - decrement to bin [" << i0 << "][" << p1 << "] - surface is " << nsurface);
                if (nsurface && nsurface != surface) break;
            }
            // find the next neighbour in Loc1
            nsurface = surface;
            while (increment(n1,bins1,open1) && surface == nsurface){
                nsurface = surfaceArrayObjects.at(n1).at(i0);
                //MSG_VERBOSE("                       - increment to bin [" << i0 << "][" << n1 << "] - surface is " << nsurface);
                if (nsurface && nsurface != surface) break;
            }
            // end th detector element
            const DetectorElementBase* element = surface->associatedDetectorElement();
            //MSG_VERBOSE("    element neighbour search for grid [" <<  p0 << " | " << i0 << " | " << n0 << "]  x [ " << p1 << " | " << i1 << " | " << n1 << " ]");
            // bools for breaking
            bool processNext1 = true; bool zero1passed = false;
            // now the previous and next bins are defined
            std::vector< const DetectorElementBase* > neighbourElements;
            for (size_t ipn1 = p1; processNext1;  increment(ipn1,bins1,open1)) {
                // zero passed for closed binning in 1
                zero1passed = zero1passed ? zero1passed : !ipn1;
                // that breaks the smaller condition
                bool processNext0 = true; bool zero0passed = false;
                for (size_t ipn0 = p0; processNext0; increment(ipn0,bins0,open0)){
                    // zero passed for closed binning in 0
                   zero0passed = zero0passed ? zero0passed : !ipn0;
                   // skip the main bin 
                   if (ipn0 == i0 && ipn1 == i1){
                       //MSG_VERBOSE("                       - skipping the bin [" << ipn0 << "][" << ipn1 << "]");
                       continue;
                   }
                   const DetectorElementBase* pnElement = surfaceArrayObjects.at(ipn1).at(ipn0) ? surfaceArrayObjects.at(ipn1).at(ipn0)->associatedDetectorElement() : nullptr;
                   //MSG_VERBOSE("                       - neighbour at bin [" << ipn0 << "][" << ipn1 << "] - element is " << pnElement);
                   if (element != pnElement && pnElement){  
                       // we can prevent double filling already here 
                       if (std::find(neighbourElements.begin(),neighbourElements.end(),pnElement) == neighbourElements.end())
                           neighbourElements.push_back(pnElement);  
                   }
                   processNext0 = (p0 < n0 ) ? ( ipn0 < n0 ) : (!zero0passed || ipn0 < n0); ;     
                }
                processNext1 = ( p1 < n1 ) ? ( ipn1 < n1 ) :  (!zero1passed || ipn1 < n1);
            }
            // standard is 8 cell connectivity, but on open borders it only knows 5 cells
            size_t neighbours         = neighbourElements.size();
            size_t expectedNeighbours = 8;
            // check the connectivity
            bool correctlyConnected = true;
            // 8-cell connectivity by default
            if (bins0 > 1 && bins1 > 1 && neighbours !=8){
                // 5-cells for the open ends 
                correctlyConnected = ( ((open0 && (p0==i0 || n0==i0)) || (open1 && (p1==i1 || n1==i1))) && neighbours ==5);
                expectedNeighbours = ( (open0 && (p0==i0 || n0==i0)) || (open1 && (p1==i1 || n1==i1)) ) ? 5 : 8;
            } else if (bins0 == 1 && neighbours != 2){
                // 2-correctlyConnected connectiviy, but 1-cell for open ends
                correctlyConnected = ( (open1 && (p1==i1 || n1==i1)) && neighbours == 1);
                expectedNeighbours =  (open1 && (p1==i1 || n1==i1)) ? 1 : 2;
            } else if (bins1 == 1 && neighbours != 2){
                // 2-cell connectiviy, but 1-cell for open ends
                correctlyConnected = ( (open0 && (p0==i0 || n0==i0)) && neighbours == 1);
                expectedNeighbours = (open0 && (p0==i0 || n0==i0)) ? 1 : 2;
            } 
            if (!correctlyConnected || !neighbours) {
                //MSG_WARNING("This element has " << neighbours << " neighbours, should be " << expectedNeighbours);
                //MSG_WARNING("    the grid for [" <<  p0 << " | " << i0 << " | " << n0 << "]  x [ " << p1 << " | " << i1 << " | " << n1 << " ]");
            }
            // register and move to the next
            element->registerNeighbours(neighbourElements);
        }
    }
    
    
}
