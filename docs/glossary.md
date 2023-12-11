# Glossary

:::{glossary}
EDM
  Event Data Model is a sta

Hole
  A hole is a state on a track where no measurement was found.

Outlier
  A measurement associated with a track, where the measurement is considered
  *far away* from the nominal parameter vector of the track on the associated
  measurement surface. This is often defined as measurement with a $\chi^2$
  larger than a certain threshold value.

SoA
  Memory layout in which individual properties of an object are stored as
  contiguous arrays in a structure, where the same index is associated with
  logically connected entries. The opposite is {term}`AoS`.

AoS
  Memory layout where a collection of objects are stored in a single array,
  where each object internally has a member for each property.
:::
