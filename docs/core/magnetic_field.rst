Magnetic field
==============

This module collects information about classes and typedefs useful for
describing different magnetic field configurations. Acts is independent of the
magnetic field implementation used. Algorithms which need magnetic field
information (e.g. :class:`Acts::RungeKuttaEngine`, :class:`Acts::AtlasStepper`,
:class:`Acts::EigenStepper`) are templated on the magnetic field. The
requirements for the magnetic field implementation are the implementation of the
following functions:

.. code-block::

    // retrieve field at given position
    Acts::Vector3D getField(const Acts::Vector3D& position) const
    // retrieve magnetic field cell at given position
    // where Cache is the specific Cache struct defined by the magnetic field implementation
    Acts::Vector3D getField(const Acts::Vector3D& position, Cache& cache) const
    // retrieve gradient of magnetic field value
    Acts::Vector3D getFieldGradient(const Acts::Vector3D& position, Acts::ActsMatrixD<3, 3>&     derivative) const
    // check whether given 3D position is inside this field cell
    bool isInside(const Acts::Vector3D& position) const

Each magnetic field implementation expects to be passed a reference to an
implementation specific cache object. This can usually be achieved through
``typename BField::Cache cache``, where ``BField`` is a template parameter. Acts
comes with the following implementations of this (implicit) interface.

Constant magnetic field implementation
--------------------------------------

Should be used to describe a constant magnetic field. The
:class:`Acts::ConstantBField` returns a given constant magnetic field value at
every point and can be set by the user either at construction or with a set
function.

Interpolated magnetic field implementation
------------------------------------------

For more complex magnetic field implementations the
:class:`Acts::InterpolatedBFieldMap` can be used. The
:class:`Acts::InterpolatedBFieldMap` internally uses a field mapper which
follows :any:`Acts::concept::AnyFieldLookup` concept. This allows users to
provide their own field mapper implementation using the
:class:`Acts::InterpolatedBFieldMap` interface. Acts provides a default field
mapper implementation :class:`Acts::InterpolatedBFieldMap::FieldMapper`, which
maps global cartesian 3D positions to magnetic field values. It uses an
underlying grid which follows the :any:`Acts::concept::AnyNDimGrid` concept which can be
a grid of any dimension and allows users to provide their own grid
implementation. Furthermore users also need to provide two functions in order to
use the :class:`Acts::InterpolatedBFieldMap::FieldMapper`:

#. A function mapping cartesian global 3D coordinates onto the grid coordinates
   of dimension N.
#. A function calculating cartesian global 3D coordinates of the magnetic field
   with the local N dimensional field and the global 3D position as an input.

A default :class:`Acts::detail::Grid` implementation is provided following the
`Acts::concept::AnyNDimGrid`, which is flexible (using template parameters) on
the dimension of the grid and the value stored in the grid. Two convenience
functions are provided to  ease the creation of an
:class:`Acts::InterpolatedBFieldMap::FieldMapper` e.g. when reading in a field
map from a file, are provided:

#. :func:`Acts::InterpolatedBFieldMap::FieldMapper<2, 2> Acts::FieldMapperRZ()`
#. :func:`Acts::InterpolatedBFieldMap::FieldMapper<3, 3> Acts::fieldMapperXYZ()`

SharedBField
------------

Wraps another ``BField`` type, which it holds as a ``std::shared_ptr<...>``. The
instance can then be copied without having to duplicate the underlying field
implementation. This is useful in the case of a large B-Field map.
