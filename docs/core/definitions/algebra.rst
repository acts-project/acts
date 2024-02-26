Algebra definitions
===================

The main algebra classes for ACTS are defined in the `Acts/Definitions/Algebra.hpp` header file.
The basic scalar type can be defined via this file and is set per default to `double`, however, if `ACTS_CUSTOM_SCALAR` is set it will be used instead.

.. code-block:: cpp

   #ifdef ACTS_CUSTOM_SCALAR
   using ActsScalar = ACTS_CUSTOM_SCALAR;
   #else 
   using ActsScalar = double;
   #endif

It is recommended within the code to deduce the Scalar type from the Event Data object, e.g.

.. code-block:: cpp

   using Scalar = Vector3::Scalar;

Currently only the `Core` package builds with `float` precision.
