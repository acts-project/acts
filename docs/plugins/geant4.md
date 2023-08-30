# Geant4 plugin

The Geant4 plugin allows to build {class}`Acts::TrackingGeometry` amnd {class}`Acts::Experimental::Detector` directly from Geant4 geometry input.
Both rely on the conversion of `G4VPhysicalVolume` into corresponding `Acts` objects.

## Surface conversion

Converting physical volumes into {class}`Acts::Surface` objects that represent sensitive detector elements, is done via the {class}`Acts::Geant4DetectorSurfaceFactory`.
This helper class allows to select volumes from the Geant4 geometry and convert them either into pairs of {class}`Acts::Geant4DetectorElement` and {class}`Acts::Surface` objects in case of sensitive elements, or simnply surfaces objects into case of passive surfaces.

The selection is hereby done by providing one or more {class}`Acts::IGeant4PhysicalVolumeSelector` objects to the surface factory.

Possible implementations of this type of conversions can be seen in the corresponding unit test `ActsUnitTestGeant4DetectorSurfaceFactory`

```c++
  // Get the box
  auto nameSelector =
      std::make_shared<Acts::Geant4PhysicalVolumeSelectors::NameSelector>(
          std::vector<std::string>{"yl"}, false);

  Acts::Geant4DetectorSurfaceFactory::Cache cache;
  Acts::Geant4DetectorSurfaceFactory::Options options;
  options.sensitiveSurfaceSelector = nameSelector;

  Acts::Geant4DetectorSurfaceFactory factory;
  factory.construct(cache, nominal, *cylinderPV, options);

  BOOST_CHECK_EQUAL(cache.sensitiveSurfaces.size(), 1u);
  BOOST_CHECK_EQUAL(cache.passiveSurfaces.size(), 0u);

  auto [ element, surface ] = cache.sensitiveSurfaces.front();
  BOOST_CHECK(surface->type() == Acts::Surface::SurfaceType::Cylinder);
```

