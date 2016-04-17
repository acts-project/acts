// STL include(s)
#include <memory>
#include <iostream>

// ACTS include(s)
#include "Tools/CylinderGeometryBuilder.h"
#include "Tools/CylinderVolumeBuilder.h"
#include "Tools/CylinderVolumeHelper.h"
#include "Tools/TrackingVolumeArrayCreator.h"
#include "Tools/LayerArrayCreator.h"
#include "Tools/PassiveLayerBuilder.h"
#include "Detector/TrackingGeometry.h"

int main()
{
  std::shared_ptr<Acts::CylinderVolumeHelper> helper(new Acts::CylinderVolumeHelper);
  helper->setLayerArrayCreator(std::unique_ptr<Acts::LayerArrayCreator>(new Acts::LayerArrayCreator));
  helper->setVolumeArrayCreator(std::unique_ptr<Acts::TrackingVolumeArrayCreator>(new Acts::TrackingVolumeArrayCreator));

  std::unique_ptr<Acts::PassiveLayerBuilder> bpLayerBuilder(new Acts::PassiveLayerBuilder);
  bpLayerBuilder->setCentralLayerRadii({21});
  bpLayerBuilder->setCentralLayerHalfLengthZ({200});
  bpLayerBuilder->setCentralLayerThickness({0.8});
  bpLayerBuilder->setCentralLayerMaterialX0({352.8});
  bpLayerBuilder->setCentralLayerMaterialL0({407.});
  bpLayerBuilder->setCentralLayerMaterialA({9.012});
  bpLayerBuilder->setCentralLayerMaterialZ({4.});
  bpLayerBuilder->setCentralLayerMaterialRho({1.848e-3});
  
  std::unique_ptr<Acts::CylinderVolumeBuilder> beamPipe(new Acts::CylinderVolumeBuilder("beamPipe"));
  beamPipe->setVolumeHelper(helper);
  beamPipe->setLayerBuilder(std::move(bpLayerBuilder));
  beamPipe->setLayerEnvelopeR(1.);
  beamPipe->setLayerEnvelopeZ(1.);

  std::unique_ptr<Acts::PassiveLayerBuilder> pixelLayerBuilder(new Acts::PassiveLayerBuilder);
  pixelLayerBuilder->setCentralLayerRadii({29,55,88,120,160,200});
  pixelLayerBuilder->setCentralLayerHalfLengthZ({5,5,5,5,5,5});
  pixelLayerBuilder->setCentralLayerThickness({0.15, 0.15, 0.15, 0.15, 0.15, 0.15});
  pixelLayerBuilder->setCentralLayerMaterialX0({95.7, 95.7, 95.7, 95.7, 95.7, 95.7});
  pixelLayerBuilder->setCentralLayerMaterialL0({465.2, 465.2, 465.2, 465.2, 465.2, 465.2});
  pixelLayerBuilder->setCentralLayerMaterialA({28.03, 28.03, 28.03, 28.03, 28.03, 28.03});
  pixelLayerBuilder->setCentralLayerMaterialZ({14, 14, 14, 14, 14, 14});
  pixelLayerBuilder->setCentralLayerMaterialRho({2.32e-3, 2.32e-3, 2.32e-3, 2.32e-3, 2.32e-3, 2.32e-3});
  pixelLayerBuilder->setPosnegLayerPositionZ({500,580,650,780});
  pixelLayerBuilder->setPosnegLayerRmin({65,65,65,65});
  pixelLayerBuilder->setPosnegLayerRmax({180,180,180,180});
  pixelLayerBuilder->setPosnegLayerThickness({0.15, 0.15, 0.15, 0.15});
  pixelLayerBuilder->setPosnegLayerMaterialX0({95.7, 95.7, 95.7, 95.7});
  pixelLayerBuilder->setPosnegLayerMaterialL0({465.2, 465.2, 465.2, 465.2});
  pixelLayerBuilder->setPosnegLayerMaterialA({28.03, 28.03, 28.03, 28.03});
  pixelLayerBuilder->setPosnegLayerMaterialZ({14, 14, 14, 14});
  pixelLayerBuilder->setPosnegLayerMaterialRho({2.32e-3, 2.32e-3, 2.32e-3, 2.32e-3});

  std::unique_ptr<Acts::CylinderVolumeBuilder> pixel(new Acts::CylinderVolumeBuilder("Pixel"));
  pixel->setVolumeHelper(helper);
  pixel->setLayerBuilder(std::move(pixelLayerBuilder));
  pixel->setLayerEnvelopeR(1.);
  pixel->setLayerEnvelopeZ(10.);
  
  std::unique_ptr<Acts::CylinderGeometryBuilder> builder(new Acts::CylinderGeometryBuilder);
  builder->setTrackingVolumeHelper(helper);
  builder->setVolumeBuilders({std::move(pixel)});
  builder->setBeamPipeBuilder(std::move(beamPipe));
  
  std::unique_ptr<Acts::TrackingGeometry> geometry = builder->trackingGeometry();
  std::cout << geometry.get() << std::endl;
  
  return 0;
}
