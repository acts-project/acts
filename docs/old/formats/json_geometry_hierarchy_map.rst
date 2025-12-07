Json Geometry Hierarchy Map
===========================

Within the code base the :class:`Acts::GeometryHierarchyMap` container is used
to map values into the geometry hierarchy. This could be e.g. track finder cuts
that are different for different parts of the detector. To simplify run-time
configuration, :class:`Acts::GeometryHierarchyMap` can be encoded/decode to/from
the following Json format. The encoded value is a Json object with two entries:
the header Json object in ``acts-geometry-hierarchy-map`` to identify the file
and value type and provide versioning for forward-compatibility, and a Json
array in ``entries`` that contains the values.

.. code-block:: json

   {
     "acts-geometry-hierarchy-map": {
       "format-version": 0,
       "value-identifier": "<user-defined-identifier>"
     },
     "entries": [
       {
         "_comment": "global default entry w/o identifier",
         "value": "..."
       },
       {
         "volume": 1,
         "layer": 2,
         "value": "..."
       },
       {
         "volume": 3,
         "value": "..."
       }
     ]
   }

Each entry is a Json object that contains the encoded content of the
:class:`Acts::GeometryIdentifier` key and the associated value. Each level
within the :class:`Acts::GeometryIdentifier` is specified by name and stored as
an integer. If a given level is not explicitly specified it is assumed to be
zero.

.. code-block:: json

   {
     "volume": "<integer>",
     "boundary": "<integer>",
     "layer": "<integer>",
     "approach": "<integer>",
     "sensitive": "<integer>",
     "value": "..."
   }

The representation of the value is not specified and depends on the specific
type.
