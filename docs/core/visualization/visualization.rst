Visualization
=============

A very lightweight layer for visualiazing Acts geometry objects and event data model is provided within the Core component.
Acts does not provide a viewer per se, but instead it was chosen to plug a visitor that can then be used for visualizing the given objects.
The visitor has to implement the `IVisualization3D` interface and can then straight forwardly used with the visualization helper structs. 
Two visualization helpers that implement industry standard 3D formats can be used from this component, 
but evidently any other visitor can be plugged in as long as it satisfies the `IVsualization` interface.

The two provided visualization visitors are:
 * `ObjVisualization3D` writing the `.obj` format, and an associated `.mtl` file for the color and material definitions
 * `PlyVisualization3D` writing the `.ply` format, which contains already the color/material information

Behind the scene
----------------

All display actions rely on the `Polyhedron` representation of Surfaces, 
i.e. each surface can be at least approximated by a list of vertices and a definition of faces connecting these vertices.
As a special feature, the `Polyhedron` can be displayed as a triangulated mesh of surfaces, i.e. each surface is divided into triangles
that build up the object to display.


Convenience helper functions
----------------------------

The visualization package also contains `static` helper functions for displaying complicated objects; these helper functions are provided by the `GeometryView3D` and `EventDataView3D` structs, that receive:
 * A visualization visitor
 * The object to be written
 * A geometry context where needed
 * Some view configuration, which is simply handled by the `ViewConfig` struct.

 The `ViewConfig` struct contains a number of parameters, such as the `visibility` flag, the RGB color definition of the object (or its contained objects),
 and some other view parameters that can be changed.

 The `Tests/UnitsTests/Core/Visualization` package contains a certain number of tests that produce output files to be visualized with any standard 3D viewer.

Example of an angular error cone:

 .. image:: ../../figures/visualization/AngularError.png
  :width: 800
  :alt: Display of an angular error. 

Example of a 2D cartesian error on a plane:

 .. image:: ../../figures/visualization/CartesianError.png
  :width: 800
  :alt: Display of an cartesian error.

Example of track parameters on a plane:

 .. image:: ../../figures/visualization/Parameters.png
  :width: 800
  :alt: Display of a track parameter object.

Example of a cylindrical layer with sensitive volumes and a surface grid:

 .. image:: ../../figures/visualization/CylinderLayer.png
  :width: 800
  :alt: Display of a track parameter object.
