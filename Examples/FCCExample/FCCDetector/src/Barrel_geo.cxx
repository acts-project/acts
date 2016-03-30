///////////////////////////////////////////////////////////////////
// Barrel_geo.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hepDetectorElement/IDetExtension.h"
#include "DD4hepDetectorElement/DetExtension.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

/**
 Constructor for a cylindrical barrel volume, possibly containing layers and the layers possibly containing modules.
*/

static Ref_t create_element(LCDD& lcdd, xml_h xml, SensitiveDetector sens)
{
    xml_det_t x_det = xml;
    string det_name = x_det.nameStr();
    //Make DetElement
    DetElement cylinderVolume(det_name, x_det.id());
    //add Extension to Detlement for the RecoGeometry
    Acts::DetExtension* detvolume = new Acts::DetExtension();
    detvolume->setShape(Acts::ShapeType::Cylinder);
    cylinderVolume.addExtension<Acts::IDetExtension>(detvolume);
    //make Volume
    DD4hep::XML::Dimension x_det_dim(x_det.dimensions());
    Tube tube_shape(x_det_dim.rmin(),x_det_dim.rmax(),x_det_dim.dz());
    Volume tube_vol(det_name,tube_shape,lcdd.air()); //air at the moment change later
    tube_vol.setVisAttributes(lcdd, x_det_dim.visStr());
    //go trough possible layers
    size_t layer_num = 0;
    for (xml_coll_t j(xml,_U(layer));j;++j)
    {
        xml_comp_t x_layer  = j;
        double l_rmin       = x_layer.inner_r();
        double l_rmax       = x_layer.outer_r();
        double l_length     = x_layer.z();
        //Create Volume and DetElement for Layer
        string layer_name  = det_name + _toString(layer_num,"layer%d");
        Volume layer_vol(layer_name,Tube(l_rmin,l_rmax,l_length), lcdd.material(x_layer.materialStr()));
        DetElement lay_det (cylinderVolume,layer_name,layer_num);
        //Visualization
        layer_vol.setVisAttributes(lcdd, x_layer.visStr());
        //go trough possible modules
        if(x_layer.hasChild(_U(module))){
            xml_comp_t x_module = x_layer.child(_U(module));
            int repeat = x_module.repeat();
            double deltaphi = 2.*M_PI/repeat;
            //slices in z
            xml_comp_t x_slice = x_layer.child(_U(slice));
            int zrepeat = x_slice.repeat();
            double dz   = x_slice.z();
            double dr   = x_slice.dr();
            size_t module_num = 0;
            //Place the Modules in z
            for (int k=-zrepeat;k<=zrepeat;k++)
            {
                double r = (l_rmax+l_rmin)*0.5;
                string zname = _toString(k,"z%d");
                if (k%2 == 0) r+=dr;
                //Place the modules in phi
                for (int i=0; i<repeat; ++i) {
                    //Creat the module volume
                    Volume mod_vol("module",Box(x_module.length(),x_module.width(),x_module.thickness()),lcdd.material(x_module.materialStr()));
                    //Visualization
                    mod_vol.setVisAttributes(lcdd, x_module.visStr());
                    double phi = deltaphi/dd4hep::rad * i;
                    string module_name = zname + _toString(i,"module%d");
                    Position trans(r*cos(phi),
                                   r*sin(phi),
                                   k*dz);
                    //Create the module DetElement
                    DetElement mod_det(lay_det,module_name,module_num);
                    //Set Sensitive Volmes sensitive
                    if (x_module.isSensitive()) {
                        mod_vol.setSensitiveDetector(sens);
                        const Segmentation segmentation(sens.readout().segmentation());
                        //add Extension for sensitive component
                        Acts::DetExtension* detSensComponent = new Acts::DetExtension();
                        detSensComponent->setSegmentation(segmentation);
                        mod_det.addExtension<Acts::IDetExtension>(detSensComponent);
                    }
                    //Place Module Box Volumes in layer
                    PlacedVolume placedmodule = layer_vol.placeVolume(mod_vol,Transform3D(RotationX(-0.5*M_PI)*RotationZ(-0.5*M_PI)*RotationX(phi-0.6*M_PI),trans));
                    placedmodule.addPhysVolID("module",module_num);
                    //assign module DetElement to the placed module volume
                    mod_det.setPlacement(placedmodule);
                    ++module_num;
                }
            }
        }
        //Place layer volume
        PlacedVolume placedLayer = tube_vol.placeVolume(layer_vol);
        placedLayer.addPhysVolID("layer",layer_num);
        //Assign layer DetElement to layer volume
        lay_det.setPlacement(placedLayer);
        ++layer_num;
    }
    //Place Volume
    Position endcap_translation(0.,0.,x_det_dim.z());
    Transform3D endcap_transform(endcap_translation);
    Volume mother_vol = lcdd.pickMotherVolume(cylinderVolume);
    PlacedVolume placedTube = mother_vol.placeVolume(tube_vol, endcap_transform);
    placedTube.addPhysVolID("system",cylinderVolume.id());
    cylinderVolume.setPlacement(placedTube);
    
    return cylinderVolume;
}

DECLARE_DETELEMENT(ACTS_Barrel, create_element)
