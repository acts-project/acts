@defgroup propagation Propagation
@brief Particle propagation

Use @ref Acts::SympyStepper or @ref Acts::EigenStepper for production propagation
in a magnetic field. The @ref Acts::AtlasStepper is a historical, hand-expanded
ATLAS Runge-Kutta-Nystroem implementation retained for performance studies. It
prioritizes speed and is not intended for production tracking. Its class
introduction records the Athena origin, numerical approximations and omitted
features.
