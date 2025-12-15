@defgroup magnetic_field Magnetic Field
@ingroup detector_descr
@ingroup propagation
@brief A description of magnetic field configurations for algorithms.

The magnetic field component of ACTS provides functionality to describe
arbitrary magnetic field configurations in an experiment. It is implemented in
a generic way and can be extended to connect to an experiment specific upstream
source of field data.

Algorithms which need magnetic field information (e.g.
@ref Acts::AtlasStepper, @ref Acts::EigenStepper) accept the magnetic
field as an explicit argument.

- The documentation of @ref Acts::MagneticFieldProvider provides a good overview
of the design and patterns associated with the magnetic field component.
- The potentially event-context dependent configuration of the magnetic field
  is documented in @ref Acts::MagneticFieldContext.

**Classes of field providers in the library:**

- Analytical magnetic fields like @ref Acts::SolenoidBField that calculate
  their respective field values on the fly.
- *Special* magnetic fields like @ref Acts::ConstantBField or @ref Acts::NullBField
- The main @ref Acts::InterpolatedBFieldMap that is most commonly used by experiments
  to provide magnetic field maps from discrete data.
