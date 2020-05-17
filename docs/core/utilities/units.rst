Unit definitions and conversions
================================

All physical quantities have both a numerical value and a unit. For the
computations we always choose a particular unit for a given physical quantity
so we only need to consider the numerical values as such. The chosen base unit
for a particular physical dimension, e.g. length, time, or energy, within this
code base is called the native unit. The base units should be chosen such that
they are internally consistent, require the least amount of explicit
conversion factors (ideally none at all), and have typical numerical values
close to unity to reduce numerical issues.

Here, the following native units are used:

- Length is expressed in mm.
- Time is expressed in [speed-of-light * time] == mm. A consequence
  of this choice is that the speed-of-light expressed in native units
  is 1.
- Angles are expressed in radian.
- Energy, mass, and momentum are all expressed in GeV (consistent with
  speed-of-light == 1).
- Electric charge is expressed in e, i.e. units of the elementary charge.
-  The magnetic field is expressed in GeV/(e*mm). The magnetic field
  connects momentum to length, e.g. in SI units the radius of a charged
  particle trajectory in a constant magnetic field is given by

  .. math::
        
     radius = - (momentum / charge) / field

  With the chosen magnetic field unit the expression above stays the
  same and no additional conversion factors are necessary.

Depending on the context a physical quantity might not be given in the native
units. In this case we need to convert to the native unit first before the value
can be used. The necessary conversion factors are defined in
:ref:`namespace_Acts__UnitConstants`. Multiplying a value with the unit constant
produces the equivalent value in the native unit, e.g.

.. code-block:: cpp

   double length = 1 * Acts::UnitConstants::m;       // length == 1000.0
   double momentum = 100 * Acts::UnitConstants::MeV; // momentum == 0.1

The conversion of a value in native units into the equivalent value in a
specific other unit is computed by dividing with the relevant unit, e.g.

.. code-block:: cpp

   double length = 10.0;                               // native units mm
   double lengthInM = length / Acts::UnitConstants::m; // == 0.01;

To further simplify the usage, physical quantities can also be expressed via
`C++ user literals <https://en.cppreference.com/w/cpp/language/user_literal>`_
define in :ref:`namespace_Acts__UnitLiterals`. This allows us to
express quantities in a concise way:

.. code-block:: cpp

   using namespace Acts::UnitLiterals;
   
   double length = 1_m;                     // == 1000.0
   double momentum = 1.25_TeV;              // == 1250.0
   double lengthInUm = length / 1_um;       // == 1000000.0
   double momentumInMeV = momentum / 1_MeV; // == 1250000.0

Since this requires a namespace import of :ref:`namespace_Acts__UnitLiterals` it
should not be used in headers since it would (accidentally) modify the namespace
wherever the header is included.

To ensure consistent computations and results the following guidelines **must**
be followed when handling physical quantities with units:

- All unqualified numerical values, i.e. without a unit, are assumed to
  be expressed in the relevant native unit, e.g. mm for lengths or GeV
  for energy/momentum.
- If a variable stores a physical quantity in a specific unit that is
  not the native unit, clearly mark this in the variable, i.e.

  .. code-block:: cpp
    
     double momentum = 100.0; // momentum is stored as native unit GeV
     double momentumInMeV = 10.0; // would be 0.01 in native units

- All input values must be given as ``numerical_value * unit_constant`` or
  equivalently using the unit literals as ``value_unit``. The resulting
  unqualified numerical value will be automatically converted to the
  native unit.
- To output an unqualified numerical value in the native units as a
  numerical value in a specific unit divide by the unit constants as
  ``numerical_value / unit_constant`` or using the unit literals as
  ``value / 1_unit``.

Examples:

.. code-block:: cpp

   #include <Acts/include/Utilities/Units.hpp>
   using namespace Acts::UnitLiterals;
   
   // define input values w/ units (via unit constants)
   double width    = 12 * Acts::UnitConstants::mm;
   double mmuon    = 105.7 * Acts::UnitConstants::MeV;
   // define input values w/ units (via unit user literals)
   double length   = 23_cm;
   double time     = 1214.2_ns;
   double angle    = 123_degree;
   double momentum = 2.5_TeV;
   double mass     = 511_keV;
   double velocity = 345_m / 1_s;
   double bfield   = 3.9_T;
   
   // convert output values (via unit constants)
   doube t_in_ns    = trackPars.time() / Acts::UnitConstants::ns;
   // convert output values (via unit user literals)
   double x_in_mm   = trackPars.position().x() / 1_mm;
   double pt_in_TeV = trackPars.momentum().pT() / 1_TeV;
