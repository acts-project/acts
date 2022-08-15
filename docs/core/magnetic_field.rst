Magnetic field
==============

.. attention::
   This section is largely **outdated** and will be replaced in the future.

This module collects information about classes and typedefs useful for
describing different magnetic field configurations. Acts is independent of the
magnetic field implementation used. Algorithms which need magnetic field
information (e.g. :class:`Acts::AtlasStepper`,
:class:`Acts::EigenStepper`) are templated on the magnetic field. The
requirements for the magnetic field implementation are shown in the example ``FieldProvider`` given below:

.. code-block:: cpp

    struct FieldProvider {
      struct Cache {
      // implementation specific
      };

      // get field for a given position, no cache
      Vector3 getField(const Vector3& pos) const;
      // get field and the gradient, no cache
      Vector3 getFieldGradient(const Vector3& pos, ActsMatrix<3,3>& deriv) const;

      // get the field for a given position, and provide the cache object
      Vector3 getField(const Vector3& position, Cache& cache) const;
      // get the field and gradient for a given position, provide cache
      Vector3 getFieldGradient(const Vector3& pos, 
                                ActsMatrix<3,3>& deriv, // mutable reference
                                Cache& cache) const;
    };

    // client code
    FieldProvider p;
    FieldProvider::Cache cache;
    auto field = p.getField({1, 2, 3}, cache); // retrieve field


Each magnetic field implementation expects to be passed a reference to an
implementation specific cache object. The cache type is provided as a nested
struct of the field provider.  It can usually be extracted using ``typename BField::Cache cache``, where ``BField`` is a template parameter. Acts comes
with the following implementations of this (implicit) interface:

- :ref:`constantbfield`
- :ref:`interpolatedbfield`
- :ref:`solenoidbfield`
- :ref:`sharedbfield`


.. _constantbfield:

Constant magnetic field
-----------------------

The simplest implementation is a constant field, which returns the same field values at every queried location. It is implemented in the :class:`Acts::ConstantBField` class.

.. doxygenclass:: Acts::ConstantBField
   :members: ConstantBField

As seen above, the class is constructed from a three-dimensional field vector, which is returned unmodified to every call to :func:`Acts::ConstantBField::getField`.

.. _interpolatedbfield:

Interpolated magnetic field
---------------------------

For more complex magnetic field implementations the
:class:`Acts::InterpolatedBFieldMap` can be used. The idea here is to calculate
an interpolated value of the magnetic field from a grid of known field values.
In 3D, this means the interpolation is done from the 8 cornerpoints of a *field
cell*. The field cell can be retrieved for any given position. Since during
typical access patterns, e.g. the propagation, subsequent steps are relatively
likely to not cross the field cell boundary, the field cell can be cached.

.. figure:: ../figures/bfield/field_cell.svg
   :width: 300

   Illustration of the field cell concept. Subsequent steps are clustered in the same field cell. The field cell only needs to be refetched when the propagation crosses into the next grid region.


The class constructor

.. doxygenfunction:: Acts::InterpolatedBFieldMap::InterpolatedBFieldMap
   :outline:

accepts a single object of type :struct:`Acts::InterpolatedBFieldMap::Config`:

.. doxygenstruct:: Acts::InterpolatedBFieldMap::Config
   :members: mapper, scale
   :outline:

The config object contains an instance of a *mapper* type, as well as a global
scale to be applied to any field values.

One implementation :struct:`Acts::InterpolatedBFieldMapper` is provided, but
since the mapper type is a template parameter, this implementation can also be
switched out. The default implementation uses :class:`Acts::detail::Grid` as
the underlying data storage. It is generic over the number of dimensions.

Most notably it exposes a type
:struct:`Acts::InterpolatedBFieldMapper::FieldCell` that corresponds to the
concept of a field cell discussed above. It also exposes a function

.. doxygenfunction:: Acts::InterpolatedBFieldMap::getFieldCell
   :outline:

that allows the retrieval of such a field cell at a given position. This function
is used by :class:`Acts::InterpolatedBFieldMap` to lookup and use field cells.
:class:`Acts::InterpolatedBFieldMap` stores the most recent field cell in
the ``Cache`` object provided by the client, and only talk to
:struct:`Acts::InterpolatedBFieldMapper` when the position leaves the current
field cell. Access to the magnetic field is done using the common interface methods

.. doxygenclass:: Acts::InterpolatedBFieldMap
   :members: getField
   :outline:

where the ``Cache`` type hides the concrete mapper used.

Helpers to construct mappers from text and root file inputs are provided:

- :func:`Acts::fieldMapperRZ`
- :func:`Acts::fieldMapperXYZ`

.. _solenoidbfield:

Analytical solenoid magnetic field
----------------------------------

Acts also provides a field provider that calculates the field vectors analytically for a solenoid field. 

.. image:: ../figures/bfield/quiver.png
   :width: 600
   :alt: Picture of a solenoid field in rz, with arrows indicating the direction of the field, and their size denoting the strength. The field is almost homogeneous in the center.

The implementation has configurable solenoid parameters:

.. doxygenstruct:: Acts::SolenoidBField::Config

.. note::

    A configuration of 

    .. code-block:: cpp

        SolenoidBField::Config cfg;
        cfg.length = 5.8_m;
        cfg.radius = (2.56 + 2.46) * 0.5 * 0.5_m;
        cfg.nCoils = 1154;
        cfg.bMagCenter = 2_T;
        SolenoidBField bField(cfg);

    roughly corresponds to the solenoid wrapping the Inner Detector in ATLAS.

Implementation
**************

The calculation uses two special functions:

- :math:`E_1(k^2)` is the complete elliptic integral of the 1st kind
- :math:`E_2(k^2)` is the complete elliptic integral of the 2nd kind

:math:`E_1(k^2)` and :math:`E_2(k^2)` are usually indicated as :math:`K(k^2)` and :math:`E(k^2)` in literature, respectively:

.. math::

  E_1(k^2) = \int_0^{\pi/2} \left( 1 - k^2 \sin^2{\theta} \right )^{-1/2} \mathop{}\!\mathrm{d}\theta

.. math::

  E_2(k^2) = \int_0^{\pi/2}\sqrt{1 - k^2 \sin^2{\theta}} \mathop{}\!\mathrm{d}\theta

:math:`k^2` is a function of the point :math:`(r, z)` and of the radius of the coil :math:`R`

.. math::

  k^2 = \frac{4Rr}{(R+r)^2 + z^2}

Using these, you can evaluate the two components :math:`B_r` and :math:`B_z` of the magnetic field:

.. math::

  B_r(r, z) = \frac{\mu_0 I}{4\pi} \frac{kz}{\sqrt{Rr^3}} \left[ \left(\frac{2-k^2}{2-2k^2}\right)E_2(k^2) - E_1(k^2) \right ]

.. math::

  B_z(r,z) = \frac{\mu_0 I}{4\pi} \frac{k}{\sqrt{Rr}} \left[ \left( \frac{(R+r)k^2-2r}{2r(1-k^2)} \right ) E_2(k^2) + E_1(k^2) \right ]

In the implementation the factor of :math:`(\mu_0\cdot I)` is defined to be a scaling factor. It is evaluated and defined as the magnetic field in the center of the coil, i.e. the scale set in :any:`Acts::SolenoidBField::Config::bMagCenter`.

.. warning::
    
   Evaluation of :math:`E_1(k^2)` and :math:`E_2(k^2)` is **slow**. The
   :class:`Acts::InterpolatedBFieldMap` easily outperforms
   :class:`Acts::SolenoidBField`. A helper :func:`Acts::solenoidFieldMapper`
   is provided that builds a map from the analytical implementation and is
   much faster to lookup.

.. _sharedbfield:

Shared magnetic field
---------------------

:class:`Acts::SharedBField` wraps another one of the magnetic field types from above.
Internally, it holds a ``std::shared_ptr<...>``, so the same field provider can be reused. This is useful in case of a larger map, for example.

.. doxygenfunction:: Acts::SharedBField::SharedBField

