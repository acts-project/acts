# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project includes
from .type_helpers import cpp_class

# python includes
from enum import Enum

""" Fundamental types """


class Type(Enum):
    SINGLE = "float"
    DOUBLE = "double"
    UINT_8 = "std::uint8_t"
    UINT_16 = "std::uint16_t"
    UINT_32 = "std::uint32_t"
    UINT_64 = "std::uint64_t"
    UINT_128 = "std::uint128_t"
    UINT_LEAST_8 = "std::uint_least8_t"
    UINT_LEAST_16 = "std::uint_least16_t"
    UINT_LEAST_32 = "std::uint_least32_t"
    UINT_LEAST_64 = "std::uint_least64_t"
    UINT_LEAST_128 = "std::uint_least128_t"

    def __str__(self) -> str:
        return self.value


""" Available linear algebra bachends """


class Algebra(Enum):
    ANY = "concepts::algebra algebra_t"
    ARRAY = "detray::array"
    EIGEN = "detray::eigen"
    FASTOR = "detray::fastor"
    SMATRIX = "detray::smatrix"
    VC_AOS = "detray::vc_aos"
    VC_SOA = "detray::vc_soa"

    def __str__(self) -> str:
        return self.value


""" Available coordinate frame types """


class Frame:
    # 2D Cartesian frame (x, y)
    CARTESIAN2D = cpp_class(specifier="detray::cartesian2D")
    # 3D Cartesian frame (x, y, z)
    CARTESIAN3D = cpp_class(specifier="detray::cartesian3D")
    # 2D concentric cylindrical frame (phi, z), no transformation
    CONCENTRIC_CYLINDRICAL2D = cpp_class(specifier="detray::concentric_cylindrical2D")
    # 2D cylindrical frame (phi, z), arbitrary placement
    CYLINDRICAL2D = cpp_class(specifier="detray::cylindrical2D")
    # 3D cylindrical frame (r, phi, z), arbitrary placement
    CYLINDRICAL3D = cpp_class(specifier="detray::cylindrical3D")
    # Linear frame (+-r, z)
    LINEAR2D = cpp_class(specifier="detray::line2D")
    # Polar frame (r, phi)
    POLAR2D = cpp_class(specifier="detray::polar2D")


""" Available geometric shape types """


class Shape:
    # ITk stereo annulus (r, phi)
    ANNULUS = cpp_class(specifier="detray::annulus2D")
    # 2D cylinder, at the origin, no rotation, single inters. sol. (r*phi, z)
    CONCENTRIC_CYLINDER = cpp_class(specifier="detray::concentric_cylinder2D")
    # 3D cubiod (x, y, z)
    CUBOID = cpp_class(specifier="detray::cuboid3D")
    # 2D cylinder (r*phi, z)
    CYLINDER2D = cpp_class(specifier="detray::cylinder2D")
    # 3D cylinder (r, phi, z)
    CYLINDER3D = cpp_class(specifier="detray::cylinder3D")
    # line, square cross sec, (+-r, z)
    DRIFT_CELL = cpp_class(specifier="detray::line_square")
    # line, circular cross sec, (+-r, z)
    STRAW_TUBE = cpp_class(specifier="detray::line_circular")
    # 2D rectangle (x, y)
    RECTANGLE = cpp_class(specifier="detray::rectangle2D")
    # 2D ring / disc (r, phi)
    RING = cpp_class(specifier="detray::ring2D")
    # 2D trapezoid (x, y)
    TRAPEZOID = cpp_class(specifier="detray::trapezoid2D")
    # Mask of any shape (to be added as 'param'), always 'inside'= true
    UNBOUNDED = cpp_class(specifier="detray::unbounded")
    # No shape
    UNMASKED = cpp_class(specifier="detray::unmasked")


""" Grid bin types """


class GridBin:
    # Single-entry bins
    SINGLE = cpp_class(specifier="detray::bins::single")
    # Satically sized bins: 9 (default)
    STATIC = cpp_class(specifier="detray::bins::static_array", param={"capacity": 9})
    # Dynamically sized bins
    DYNAMIC = cpp_class(specifier="detray::bins::dynamic_array")


""" Grid serializer types """


class GridSerializer:
    # Single-entry bins
    SIMPLE = cpp_class(specifier="detray::simple_serializer")


MATERIAL_MAP_SPECIFIER = "detray::material_map"

""" Available material types """


class Material:
    # Slab of material of given thickness
    SLAB = cpp_class(specifier="detray::material_slab")
    # Material with round cross sec, of given r
    ROD = cpp_class(specifier="detray::material_rod")
    # Raw material type, used e.g. for homogeneous volume material
    RAW = cpp_class(specifier="detray::material")
    # Surface material map, annulus shape 2D
    ANNULUS_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.ANNULUS, "frame": Frame.POLAR2D},
    )
    # Surface material map, concentric cyl.
    CONCENTIRC_CYLINDER_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={
            "shape": Shape.CONCENTRIC_CYLINDER,
            "frame": Frame.CONCENTRIC_CYLINDRICAL2D,
        },
    )
    # Surface material map, cylindrical 2D
    CYLINDER_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.CYLINDER2D, "frame": Frame.CYLINDRICAL2D},
    )
    # Volume material map, cylindrical 3D
    CYLINDER_MAP3D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.CYLINDER3D, "frame": Frame.CYLINDRICAL3D},
    )
    # Volume material map, cuboid 3D
    CUBOID_MAP3D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.CUBOID, "frame": Frame.CARTESIAN3D},
    )
    # Surface material map, rectangular 2D
    RECTANGLE_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.RECTANGLE, "frame": Frame.CARTESIAN2D},
    )
    # Surface material map, disc 2D
    DISC_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.RING, "frame": Frame.POLAR2D},
    )
    # Surface material map, trapezoidal 2D
    TRAPEZOID_MAP2D = cpp_class(
        specifier=MATERIAL_MAP_SPECIFIER,
        param={"shape": Shape.TRAPEZOID, "frame": Frame.CARTESIAN2D},
    )


SPATIAL_GRID_SPECIFIER = "detray::spatial_grid"

""" Available surface/volume acceleration structure types """


class Accelerator:
    # Test all registered surfaces
    BRUTE_FORCE = cpp_class(specifier="detray::brute_force")
    # Surface grid, cylindrical 2D (barrel)
    CONCENTRIC_CYLINDER_GRID2D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.CONCENTRIC_CYLINDER,
            "frame": Frame.CONCENTRIC_CYLINDRICAL2D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
    # Surface grid, cylindrical 2D (barrel)
    CYLINDER_GRID2D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.CYLINDER2D,
            "frame": Frame.CYLINDRICAL2D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
    # Surface grid, cylindrical 3D (barrel)
    CYLINDER_GRID3D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.CYLINDER3D,
            "frame": Frame.CYLINDRICAL3D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
    # Surface grid, disc (endcap)
    DISC_GRID2D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.RING,
            "frame": Frame.POLAR2D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
    # Surface grid, reactangular (telescope)
    RECTANGLE_GRID2D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.RECTANGLE,
            "frame": Frame.CARTESIAN2D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
    # Surface grid, cuboid 3D (telescope)
    CUBOID_GRID3D = cpp_class(
        specifier=SPATIAL_GRID_SPECIFIER,
        param={
            "shape": Shape.CUBOID,
            "frame": Frame.CARTESIAN3D,
            "bin": GridBin.DYNAMIC,
            "serializer": GridSerializer.SIMPLE,
        },
    )
