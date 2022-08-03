.. _dd4hep_old:

DD4hep plugin (old integration model)
=============

.. warning::
   The DD4hep integration mechanism has been removed in `v20.0.0
   <https://github.com/acts-project/acts/releases/tag/v20.0.0>`_. Starting from
   this version, the DD4hep plugin requires a minimum DD4hep version of
   ``1.21``, which includes facilities for a new integration model.

   See :ref:`DD4hep plugin` for documentation on the new integration model.

The DD4hepPlugin allows building of a :class:`Acts::TrackingGeometry` from
`DD4hep`_ input. DD4hep uses `ROOT`_ TGeo as the underlying geometry model.

.. _DD4hep: https://dd4hep.web.cern.ch/dd4hep/
.. _ROOT: https://root.cern.ch

General
-------

The basic input for building the detector are a detector description in the
XML file format and a corresponding detector constructor written in C++. These
two components have to be changed accordingly. Detector constructors use these
XML files as an input and construct a detector in the DD4hep geometry format.
The whole detector is segmented into different detector parts, e.g. barrels
and endcaps which describe different sub-detectors. Since these sub-detectors
are built differently, they need different detector constructors. In this way
one detector description in XML can use various detector constructors, needed
for the different kinds of detector parts.

In the detector description model of DD4hep any detector is a tree of instances
of the so-called ``DetElement`` class. This ``DetElement`` class provides all
needed detector information, e.g. readout segmentation, geometrical information,
environmental conditions. This tree is parallel to the volume tree, which
provides the ``TGeoVolumes`` and their placements. The relation between these
two trees is one-directional, i.e. every volume can be accessed via its
corresponding ``DetElement``, but not vice versa. Not every supporting material
will be declared as a detector element, hence, the geometrical tree can have a
deeper hierarchy structure. In order to access both, detector specific and
geometrical information, the conversion to the tracking geometry navigates
through the detector tree. The ``DetElement`` can also be extended, to add
specific features or to access information. This extension mechanism is used
during the translation process.

Acts extension
--------------

.. warning::
   The DD4hep integration mechanism has been removed in `v20.0.0
   <https://github.com/acts-project/acts/releases/tag/v20.0.0>`_. Starting from
   this version, the DD4hep plugin requires a minimum DD4hep version of
   ``1.21``, which includes facilities for a new integration model.

DD4hep provides a special extension mechanism for the ``DetElement`` which
allows to add custom features. In Acts this functionality is used for the
conversion from ``DD4hep`` into Acts. The extensions are used to indicate
certain volumes, e.g. if a ``DetElement`` is the beam pipe or if a
``DetElement`` is a layer carrying the sensitive modules. In addition the
extensions are used in order to distinguish if a sub detector is a barrel
(described as a cylinder volume in Acts) or an endcap (which is described as a
disc volume in Acts). Furthermore the extensions are used to hand over specific
information needed for tracking, e.g. paramters for material mapping.

DD4hepDetectorElement
---------------------

In Acts the surfaces describing the sensitive modules of a detector are directly
linked to these of the initial geometry input. In the case of DD4hep the
:class:`Acts::DD4hepDetectorElement` was introduced which is the direct link of
Acts to DD4hep. In the case for tracking relevant parameters in the DD4hep
geometry description are changed (e.g. alignment) it will be automatically
changed in Acts.

Build
-----

The DD4hepPlugin is only build on demand. The DD4hepPlugin depends on the
TGeoPlugin therefore both plugins need to be installed. During the cmake
configuration the flags ``ACTS_BUILD_PLUGIN_DD4HEP=on`` and
``ACTS_BUILD_PLUGIN_TGEO=on`` need to be set. In addition ROOT and DD4hep
installations need to be available to cmake.

Prerequisites
-------------

To guarantee a working translation from DD4hep input to Acts geometry the
following conditions need to be met:

- The detector needs to have a barrel-endcap structure: Every hierarchy of
  subdetectors (e.g. PixelDetector, StripDetector,...) needs to be decomposed
  into
  
  #. {barrel}
  #. {barrel + 2 endcaps}
  #. {2 endcaps} - in case there is no barrel at this stage (e.g. forward end caps)

- If a hierachy is not only a single barrel but is decomposed of a barrel
  and its corresponding endcaps they need to be grouped together in an
  assembly using the ``DD4hep_SubdetectorAssembly`` constructor which is
  provided by DD4hep. Example of usage in xml file (where Barrel0, nEndCap0
  and pEndCap0 are sub detectors defined in the file ``PixelTracker.xml``):
  
  .. code-block:: xml
  
     <include ref="PixelTracker.xml"/>
     <detectors>
       <detector id="1" name="PixelTracker"type="DD4hep_SubdetectorAssembly"vis="BlueVisTrans">
         <shape name="PixelEnvelope" type="Tube" rmin="Env0_rmin"rmax="Env0_rmax"dz="Env0_dz" material="Air"/>
           <composite name="Barrel0"/>
           <composite name="nEndCap0"/>
           <composite name="pEndCap0"/>
       </detector>
     </detectors>

  If a user wants to create his/her own constructor to group these
  volumes together the type needs to be set to "compound".

- Since the translation walks trough the ``DetElement`` tree the following
  objects need to be declared as a DD4hep ``DetElement``:
 
  - The subvolumes e.g. barrel, endcap, beampipe (they are usually build with
    different DD4hep constructors and are therefore DD4hep ``DetElement``'s
    per default).
  - Layers when containing sensitive material and/or the layer should
    carry material (which will be mapped on the layer if indicated), or
    the layer is sensitive itself.
  
    .. note::
    
       the layer does not need to be a direct child of the volume (barrel or
       endcap),it an be nested in substructures

  - Sensitive detector modules
    
    .. note::
      
       The sensitive detector modules need to be placed in a layer however
       it can be nested in substructures (can be a component of a modules)
       i.e. it does not need to be a direct child of the layer

- The Tracking geometry needs to be built from bottom to top to ensure
  navigation. Therefore the different hierarchies need to be sorted ascending.
  Per default the sub detectors are sorted by the id of their ``DetElement``.
  In case another sorting needs to be applied, the users can provide their own
  function.

- The :class:`Acts::ActsExtension`'s need to be used during the detector
  construction indicating if a ``DetElement``
  
  - is a barrel
  - is an endcap
  - is the beampipe
  - is a layer

There are two modes building the layers around the sensitive detector modules:

- The ``DetElement`` containing the sensitive modules have a geometrical
  shape.
  
  The boundaries of the layers in Acts are taken directly from the given shape.

- The ``DetElement`` containing the sensitive modules have no specific shape
  (assembly).
  
  The boundaries of the layers are calculated automatically by adding a
  tolerance to the geometric extension of the contained surfaces. The
  tolerances in r and z need to be set for every ``DetElement`` representing
  layer using envelopeR and envelopeZ in :class:`Acts::ActsExtension`.

The volumes are automatically build around the layers:

- The boundaries for the volumes are calculated automatically by adding a
  tolerance to the geometric extension of the contained layers. The
  tolerance parameters ``layerEnvelopeR`` and ``layerEnvelopeZ`` need to be
  set in the :func:`Acts::convertDD4hepDetector()` function.

Furthermore parameters can be handed over for material mapping or the axes
orientation of modules.

Summing up the ``DetElement`` tree in DD4hep should have the following
structure:

.. image:: /figures/DD4hepPlugin_DetElementStructure.jpg

It is also possible to translate a very simple detector geometry, which just
consists of cylindrical (for a barrel) or disc (for endcaps) layers which either
have material, or, are declared sensitive in dd4hep themselves without
containing any detector modules.

Usage
-----

To receive the :class:`Acts::TrackingGeometry` the user should use the global
function :func:`Acts::convertDD4hepDetector()`, where he/she needs to hand over
the DD4hep world ``DetElement``. For a valid translation the user needs to make
sure, that all prerequisites described above are met and that the right
:class:`Acts::ActsExtension`'s are added during the DD4hep construction.
