"""GeometryIdentifier bit layout (see ``Acts::GeometryIdentifier`` in C++)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GeometryIdentifierMasks:
    """Masks and shifts match ``Core/include/Acts/Geometry/GeometryIdentifier.hpp``."""

    volume_mask: int = 0xFF00000000000000
    boundary_mask: int = 0x00FF000000000000
    layer_mask: int = 0x0000FFF000000000
    approach_mask: int = 0x0000000FF0000000
    sensitive_mask: int = 0x000000000FFFFF00
    extra_mask: int = 0x00000000000000FF

    @property
    def volume_shift(self) -> int:
        return (self.volume_mask & -self.volume_mask).bit_length() - 1

    @property
    def boundary_shift(self) -> int:
        return (self.boundary_mask & -self.boundary_mask).bit_length() - 1

    @property
    def layer_shift(self) -> int:
        return (self.layer_mask & -self.layer_mask).bit_length() - 1

    @property
    def approach_shift(self) -> int:
        return (self.approach_mask & -self.approach_mask).bit_length() - 1

    @property
    def sensitive_shift(self) -> int:
        return (self.sensitive_mask & -self.sensitive_mask).bit_length() - 1

    @property
    def extra_shift(self) -> int:
        return (self.extra_mask & -self.extra_mask).bit_length() - 1


DEFAULT_MASKS = GeometryIdentifierMasks()
