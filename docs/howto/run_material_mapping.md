# ACTS Material Mapping Tutorial
This page will explain how to perform the material mapping with the ACTS Examples. For this example we will use the Open Data Detector (ODD) the last paragraph will explain what need to be change if you want to perform the material mapping with another detector.

## Prerequisites
First you will need to build the ACTS with the Examples, Geant4 and the json plugin (`ACTS_BUILD_EXAMPLES`, `ACTS_BUILD_EXAMPLES_GEANT4` and `ACTS_BUILD_JSON_PLUGIN`), please refer to the general how-to ACTS guide. Depending on the type of detector you want to map you will need to use some additional package, in our case `ACTS_BUILD_EXAMPLES_DD4HEP` and `ACTS_BUILD_TGEO_PLUGIN` are needed.

Once the Acts has been built we can start the mapping. The mapping is divided in two aspects : the surface mapping in which the material is mapped onto the closest surfaces (following the propagation direction) and the volume mapping in which the material is mapped onto a 3D (or 2D) grid associated to a volume. The first step is to select which surfaces and volumes we will want to map material onto, this is done by association a ProtoSurfaceMaterial (or a ProtoVolumeMaterial) to the surfaces (or volumes) of interest. In the case of the ODD and some other DD4hep detectors this is done at the building step. For other detector or if one want to be able to control precisely which layer will be mapped on and with which binning an additional step is required.

## Mapping and configuration

First we will need to extract the list of all the surfaces and volumes in our detector, to do so we will use the GeometryExample :

```
./../build/bin/ActsExampleGeometryDD4hep -n1 -j1 --mat-output-file geometry-map  --dd4hep-input ../Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml --output-json true --mat-output-allmaterial true --mat-output-sensitive false
```

This Algorithm is useful to obtain a visualisation of your detector using the different type of output available (`output-obj` give .obj with 3D representation of you different subdetector for example), here we use `output-json` to obtain a map of all the surfaces and volumes in the detector with a ProtoSurfaceMaterial (or a ProtoVolumeMaterial), `mat-output-allmaterial` ensure that a ProtoSurfaceMaterial (or a ProtoVolumeMaterial) is associated to all the surfaces (or volumes) enforcing all of them to be written.
Four type of surfaces exist :  
*  boundaries which represent the boundaries of the different volumes
*  approaches which to respond to the entrance and exist of the detector layer
*  representing which correspond to the surface representation of a layer (often positioned at the middle of the 2 approaches)
*  sensitives which correspond to the active part of the detector (sensors)

By default, all the surfaces will be written but one can turn a specific type off (for example the sensitive) by using the appropriate option : `mat-output-XXX false`

The Json file can now be edited to select which surfaces and volumes you want to have material mapped on. The Json file contain a list of volumes, in those volumes you can find boundaries and layers, the layers then contains approaches, representing and sensitives surfaces. Information of the surface such as their type, range, id and position are available. To add one surface to the material mapping one simply need to switch the `mapMaterial` variable to true, the binning can then be change by changing the number associated to `bin0` and `bin1` the type of bin can also be change. For the volume the same method can be applied except that up to 3 bins can be associated.

When mapping onto a surface, the material inside volumes with material (or ProtoMaterial) will be ignored, you should thus avoid mapping material onto surface within material volumes. When mapping onto a volume, only the material within that volume will be used. If you have a large gap between you last material surface and your volume you might then want to also map material onto the boundary of the material volume.


This configuration can be cumbersome to do by hand especially when trying to map on sensitives surfaces. To simplify this task two python script are available in Examples/scripts/MaterialMaping :

*  writeMapConfig.py
*  configureMap.py

The first one take as an input the surfaces map previously generated and will return a json config file. In this file you can see all the different type of surfaces associated to each volume. You can then change the binning associated to a surface type, when the second script is called the resulting map will have the binning applied to all the surface of the corresponding type and `"mapMaterial"` will be change to true. Attention : the configureMap.py will modify the surfaces map used in input you might want to make a copy of it beforehand.


```
python3 ../acts-framework/scripts/MaterialMaping/writeMapConfig.py geometry-map.json config-map.json
```

Then edit the config-map.json file

```
../acts-framework/scripts/MaterialMaping/configureMap.py geometry-map.json config-map.json
```

## Geantino scan

The next step is to do a geantino scan of our detector. For this we will use the GeantinoRecording algorithm :

```
./../build/bin/ActsSimGeantinoRecording -j1 --dd4hep-input ../Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml --output-root true -n10000
```

The result of the geantino scan will be a root file containing Material Tracks, those contain the direction and production vertex of the geantino, the total material accumulated and all the interaction point in the detector.

## Material Mapping

With the surfaces map and the material track we can finally do the material mapping using the MaterialMapping algorithm :

```
./../build/bin/ActsExampleMaterialMappingDD4hep -j1 --input-root true --input-files geant-material-tracks.root --mat-input-type file --mat-input-file geometry-map.json --output-root true --output-json true --mat-mapping-collection material-tracks --mat-output-file material-maps --mat-mapping-surfaces true --mat-mapping-volumes true --mat-mapping-volume-stepsize 1 --dd4hep-input ../Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml
```

As an output you will obtain the material map as a root and json file and a new Material Tracks collection in a root file. This new collection add to each material interaction the associated surface during the mapping. This can be used for the control plots.
Depending on what you want to do there are three option you can change :
*  `mat-mapping-surfaces` : determine if material is mapped onto surfaces
*  `mat-mapping-volumes` : determine if material is mapped onto volumes
*  `mat-mapping-volume-stepsize` : Determine the step size used in the sampling of the volume this should be small compared to you bin size.

You can map onto surfaces and volumes separately (for example if you want to optimise one then the other). In that case after mapping one of those you will need to use the resulting json material map as an input to the `mat-input-file`.

## Material Validation

Now that the map has been written, you will want to validate it. First you can use the MaterialValidation example. This will perform propagation throughout the detector once it has been decorated with the material map. It will then output material tracks with the same format as the one obtain with the Geantino.

```
./../build/bin/ActsExampleMaterialValidationDD4hep -n 1000 --mat-input-type file --mat-input-file material-maps.json --output-root true --mat-output-file val-mat-map --dd4hep-input ../Examples/Detectors/DD4hepDetector/compact/OpenDataDetector/OpenDataDetector.xml
```

To do the validation, five root macro are available in scripts/MaterialMaping :

*  Mat_map.C : General comparison at the track level and 2D map of the detector.
*  Mat_map_surface_plot.C : for each mapped surface show the position of the material.
*  Mat_map_surface_plot_ratio.C : material ration between the truth and the validation for each surface.
*  Mat_map_surface_plot_dist.C : position of the geantino interaction with respect to the surface they are mapped on.
*  Mat_map_surface_plot_1D.C : 1D distribution of the material in each surface.

```
mkdir Validation

root -l -b ../acts-framework/scripts/MaterialMaping/Mat_map.C'("propagation-material.root","material-maps_tracks.root","Validation")'
.q

mkdir Surfaces
cd Surfaces
mkdir prop_plot
mkdir map_plot
mkdir ratio_plot
mkdir dist_plot
mkdir 1D_plot
cd ..

root -l -b ../acts-framework/scripts/MaterialMaping/Mat_map_surface_plot_ratio.C'("propagation-material.root","material-maps_tracks.root","geometry-map.json",100000,"Surfaces/ratio_plot","Surfaces/prop_plot","Surfaces/map_plot")'
.q
root -l -b ../acts-framework/scripts/MaterialMaping/Mat_map_surface_plot_dist.C'("material-maps_tracks.root","geometry-map.json",-1,"Surfaces/dist_plot")'
.q
root -l -b ../acts-framework/scripts/MaterialMaping/Mat_map_surface_plot_1D.C'("material-maps_tracks.root","geometry-map.json",100000,"Surfaces/1D_plot")'
.q
```

Using the validation plot you can then adapt the binning and the mapped surface to improve the mapping.
Depending on your root version those macros might not work. They have been tested with version 6.18.04 so you can always revert to that version in case of problem.

## Using a different detector

If you want to use a different detector type of detector you will first need to ensure that the relevant packages were added during the compilation. After that if your detector is a DD4hep detector you will just need to replace the path given to the '--dd4hep-input' option. In case it is another type of detector implementation you can replace DD4hep in the name of the algorithm by what correspond to your detector implementation. For more information on how to include your detector in that case you can refer to the documentation of the algorithm using the `-h` option.
