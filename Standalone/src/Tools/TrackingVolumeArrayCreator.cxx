///////////////////////////////////////////////////////////////////
// TrackingVolumeArrayCreator.cxx, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Tools/TrackingVolumeArrayCreator.h"
#include "ACTS/GeometryUtils/BinUtility.h"
#include "ACTS/GeometryUtils/BinnedArray1D.h"
#include "ACTS/GeometryUtils/GeometryObjectSorter.h"
#include "ACTS/Detector/TrackingVolume.h"
#include "ACTS/Volumes/VolumeBounds.h"
// Core module
#include "ACTS/Utilities/AlgebraDefinitions.h"


std::shared_ptr<const Acts::TrackingVolumeArray> Acts::TrackingVolumeArrayCreator::trackingVolumeArray(const TrackingVolumeVector& tVolumes, BinningValue bValue) const 
{
    // MSG_VERBOSE("Create VolumeArray of "<< tVolumes.size() << " TrackingVolumes with binning in : " << binningValueNames[bValue] );
    // let's copy and sort 
    TrackingVolumeVector volumes(tVolumes);
    // sort it accordingly to the binning value
    GeometryObjectSorterT< std::shared_ptr<const TrackingVolume> > volumeSorter(bValue);
    std::sort(volumes.begin(),volumes.end(),volumeSorter);
    
    // prepare what we need :
    // (1) arbitrary binning for volumes is fast enough
    std::vector<float> boundaries;
    boundaries.reserve(tVolumes.size()+1);
    // (2) the vector needed for the BinnedArray
    std::vector<TrackingVolumeOrderPosition> tVolumesOrdered;
        
    // let's loop over the (sorted) volumes
    for (auto& tVolume : volumes ){
        // get the binning position
        Vector3D binningPosition = tVolume->binningPosition(bValue);
        double binningBorder = tVolume->volumeBounds().binningBorder(bValue);
        // get the center value according to the bin
        double value = tVolume->binningPositionValue(bValue);
        // for the first one take low and high boundary
        if (!boundaries.size()) boundaries.push_back(value-binningBorder);
        // always take the high boundary
        boundaries.push_back(value+binningBorder);
        // record the volume to be ordered 
        tVolumesOrdered.push_back(TrackingVolumeOrderPosition(tVolume,binningPosition));
    }
    
    // now create teh bin utility
    BinUtility* binUtility = new BinUtility(boundaries, open, bValue);

    // and return the newly created binned array
    return std::shared_ptr<const TrackingVolumeArray>(new BinnedArray1D< std::shared_ptr<const TrackingVolume> >(tVolumesOrdered,binUtility));
}
