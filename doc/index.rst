A Common Tracking Software
==========================


This project contains an experiment-independent set of track reconstruction
tools. The main philosophy is to provide high-level track reconstruction modules
that can be used for any tracking detector. The description of the tracking
detector's geometry is optimized for efficient navigation and fast
extrapolation of tracks. Converters for several common geometry description
languages exist. Having a highly performant, yet largely customizable
implementation of track reconstruction algorithms was a primary objective for
the design of this toolset. Additionally, the applicability to real-life HEP
experiments plays major role in the development process. Apart from algorithmic
code, this project also provides an event data model for the description of
track parameters and measurements.

Key features of this project include:

*   tracking geometry description which can be constructed from TGeo, DD4Hep, or GDML input,
*   simple and efficient event data model,
*   performant and highly flexible algorithms for track propagation and fitting,
*   basic seed finding algorithms.


.. toctree::
    :maxdepth: 2

    modules/modules
    plugins/plugins
    api/api
    grid_axis
    logging
