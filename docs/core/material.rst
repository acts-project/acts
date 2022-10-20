.. _material_core:

Material
========

.. note::

   This was originally written for the materials plugin that has been merged
   into the material module

The material module allow to map material from a detailed full detector geometry
onto the simplified Acts geometry. The material is mapped onto layers of the
tracking geometry which are marked to carry support material. The marking is
done during the geometry building process. The material can be mapped onto
either, the inner, the outer boundary surface or the central (representing)
:class:`Acts::Surface` of the :class:`Acts::Layer`. The :class:`Acts::Material`
is described on a two dimensional grid for each layer
(:class:`Acts::BinnedSurfaceMaterial`). The user defines the granularity of the
grid during the geometry building process.

.. note::

   The :doc:`/plugins/dd4hep` offers the possibility to mark layers which should
   carry material and to determine the grid granularity, using the class
   :class:`Acts::ActsExtension`.

Following the Acts philosophy the material mapping is agnostic to any file
format and software used to create or store the material maps. The material
should be stored in instances of the class :class:`Acts::MaterialTrack`. This
material track record represents a track starting from a certain position, in a
certain direction, containing all material along this track. The material along
the material track record is stored as a container of
:class:`Acts::MaterialStep` instances. Each material step contains the material
and its thickness at a certain position.

The material mapping process can be split into two subprocesses:

- Material assignment
- Material averaging

Material assignment
-------------------

During the material assignment process the decision onto which layer each
material step will be assigned is done. To assign a :class:`Acts::MaterialTrack`
the function :func:`Acts::MaterialMapping::mapMaterial()` should be used. This
function extrapolates through the tracking detector, with the start position and
direction given by the material track record and collects all layers marked to
carry material. Then it loops through all material steps of the material track
record and assigns the material of each step to the closest layer:

.. figure:: /figures/MaterialAssignment.jpeg

   Example of material assignment onto the inner boundary surface of the
   layers. The green points are assigned to the current inner layer, the red to
   the next inner layer."

Material averaging
------------------

During the material mapping the user can decide to average the material whenever
they prefer by using the function
:func:`Acts::MaterialMapping::averageLayerMaterial()`. In the end when all
material track records have been mapped one should use the function
:func:`Acts::MaterialMapping::finalizeLayerMaterial()` in order to finalize the
process.

The full material mapping process should be done in the framework of the user.

Possible workflow:

- Create material map(s) of full detector geometry using :class:`Acts::MaterialTrack`.
- Read in material map(s) and go through all collected material track records.
  Use :func:`Acts::MaterialMapping::mapMaterial()` for each material track
  record.
- Use :func:`Acts::MaterialMapping::averageLayerMaterial()` once per run.
- In the end of the process use
  :func:`Acts::MaterialMapping::finalizeLayerMaterial()` which assigns the
  final material to the layers.
