Acts Common Tracking Software
=============================

Acts is an experiment-independent toolkit for (charged) particle track
reconstruction in (high energy) physics experiments implemented in modern C++.

The Acts project provides high-level track reconstruction modules that can be
used for any tracking detector. The tracking detector geometry description is
optimized for efficient navigation and fast extrapolation of tracks. Converters
for several common geometry description packages are available. In addition to
the algorithmic code, this project also provides an event data model for the
description of track parameters and measurements.

Key features:

* A tracking geometry description which can be constructed manually or from
  TGeo and DD4hep input.
* Simple event data model.
* Implementations of common algorithms for track propagation and fitting.
* Implementations of basic seed finding algorithms.
* Implementations of common vertexing algorithms.

.. toctree::
   :maxdepth: 3

   getting_started
   tracking
   acts_project
   core/core
   fatras/fatras
   plugins/plugins
   examples/examples
   howto/howto
   contribution/contribution

   versioning
   formats/formats
   api/api
   codeguide
   authors
   license
