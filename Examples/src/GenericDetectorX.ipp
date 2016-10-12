//-------------------------------------------------------------------------------------
// Pixel volume
//-------------------------------------------------------------------------------------

// configure pixel layer builder
GenericLayerBuilder::Config layerBuilderConfigPixel;
layerBuilderConfigPixel.layerCreator        = layerCreator;
layerBuilderConfigPixel.layerIdentification = "Pixel";

// central layers:

layerBuilderConfigPixel.centralLayerBinMultipliers = {1, 1};

// fill necessary vectors for configuration
//-------------------------------------------------------------------------------------

std::vector<double> centralLayerStaggerZ;
std::vector<double> centralLayerOverlapZ;

centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
centralLayerStaggerZ.push_back(0.5);
centralLayerOverlapZ.push_back(0.625);
layerBuilderConfigPixel.centralLayerRadii = {29, 88, 160, 230, 320, 440};
layerBuilderConfigPixel.centralLayerEnvelopes = {{0.5, 2.0}, {0.5, 2.0}, {0.5, 2.0}, {0.5, 2.0}, {0.5, 2.0}, {0.5, 2.0}};
layerBuilderConfigPixel.centralLayerMaterialConcentration = {1, 1, 1, 1, 1, 1};
layerBuilderConfigPixel.centralModuleBinningSchema = {{13,25}, {39,25}, {72,25}, {104,25}, {145,25}, {199,25}};
layerBuilderConfigPixel.centralModuleTiltPhi = {0.18, 0.18, 0.18, 0.18, 0.18, 0.18};
layerBuilderConfigPixel.centralModuleHalfX = {8.3, 8.3, 8.3, 8.3, 8.3, 8.3};
layerBuilderConfigPixel.centralModuleHalfY = {16.3, 16.3, 16.3, 16.3, 16.3, 16.3};
layerBuilderConfigPixel.centralModuleThickness = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
// mPositions
std::vector<std::vector<Vector3D>> centralModulePositions;
for (size_t cl = 0; cl < layerBuilderConfigPixel.centralLayerRadii.size(); ++cl) {
  // call the helper function
  centralModulePositions.push_back(
      modulePositionsCylinder(layerBuilderConfigPixel.centralLayerRadii[cl],
                              centralLayerStaggerZ[cl],  // staggering in z
                              layerBuilderConfigPixel.centralModuleHalfY[cl],
                              centralLayerOverlapZ[cl],  // 2 mm module overlap
                              layerBuilderConfigPixel.centralModuleBinningSchema[cl]));
}
layerBuilderConfigPixel.centralModulePositions = centralModulePositions;

// central layers:

// define the builder
auto layerBuilderPixel = std::make_shared<GenericLayerBuilder>(
   layerBuilderConfigPixel, getDefaultLogger("PixelVolumeBuilder", layerLLevel));
// build the Pixel volume
CylinderVolumeBuilder::Config volumeBuilderConfigPixel;
volumeBuilderConfigPixel.trackingVolumeHelper = cylinderVolumeHelper;
volumeBuilderConfigPixel.volumeName           = "Pixel";
volumeBuilderConfigPixel.buildToRadiusZero    = false;
volumeBuilderConfigPixel.layerEnvelopeR       = {1. * Acts::units::_mm, 5. * Acts::units::_mm};
volumeBuilderConfigPixel.layerBuilder         = layerBuilderPixel;
volumeBuilderConfigPixel.volumeSignature      = 0;
auto volumeBuilderPixel        = std::make_shared<CylinderVolumeBuilder>(
    volumeBuilderConfigPixel, getDefaultLogger("PixelVolumeBuilder", volumeLLevel));

// add to the list of builders
volumeBuilders.push_back(volumeBuilderPixel);