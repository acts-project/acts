# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

"""Schema of a homogeneous material file"""

homogeneous_material_schema = {
    "type": "object",
    "properties": {
        "header": {
            "type": "object",
            "description": "File header section with important IO metadata",
            "properties": {
                "common": {
                    "type": "object",
                    "description": "General metadata",
                    "properties": {
                        "version": {
                            "type": "string",
                            "description": "detray version this file was created with",
                        },
                        "detector": {
                            "type": "string",
                            "description": "Name of the detector",
                        },
                        "date": {
                            "type": "string",
                            "description": "Time stamp (file creation)",
                        },
                        "tag": {
                            "type": "string",
                            "description": "Type of data",
                            "enum": ["homogeneous_material"],
                        },
                    },
                    "required": ["version", "detector", "tag"],
                },
                "slab_count": {
                    "type": "integer",
                    "description": "Number of material slabs in the detector",
                    "minimum": 0,
                },
                "rod_count": {
                    "type": "integer",
                    "description": "Number of material rods in the detector",
                    "minimum": 0,
                },
            },
            "required": ["common"],
        },
        "data": {
            "type": "object",
            "description": "Homogeneous material data (volume-by-volume)",
            "properties": {
                "volumes": {
                    "type": "array",
                    "minitems": 1,
                    "items": {
                        "type": "object",
                        "description": "Material data of a single volume",
                        "properties": {
                            "volume_link": {
                                "type": "integer",
                                "description": "The volume index",
                            },
                            "material_slabs": {
                                "type": "array",
                                "minitems": 1,
                                "items": {
                                    "type": "object",
                                    "description": "The material slab",
                                    "properties": {
                                        "type": {
                                            "type": "integer",
                                            "description": "Material type ID",
                                            "enum": [5],
                                        },
                                        "surface_idx": {
                                            "type": "integer",
                                            "description": "Volume-local index of the surface the slab belongs to",
                                            "minimum": 0,
                                        },
                                        "thickness": {
                                            "type": "number",
                                            "description": "Thickness of the slab",
                                            "minimum": 0,
                                        },
                                        "index_in_coll": {
                                            "type": "integer",
                                            "description": "Volume-local sorting index",
                                            "minimum": 0,
                                        },
                                        "material": {
                                            "type": "object",
                                            "properties": {
                                                "params": {
                                                    "type": "array",
                                                    "description": "Material parameters",
                                                    "minitems": 7,
                                                    "maxItems": 7,
                                                    "items": {"type": "number"},
                                                }
                                            },
                                            "required": ["params"],
                                        },
                                    },
                                    "required": [
                                        "type",
                                        "surface_idx",
                                        "thickness",
                                        "material",
                                    ],
                                },
                            },
                            "material_rods": {
                                "type": "array",
                                "minitems": 1,
                                "items": {
                                    "type": "object",
                                    "description": "The material rods",
                                    "properties": {
                                        "type": {
                                            "type": "integer",
                                            "description": "Material type ID",
                                            "enum": [6],
                                        },
                                        "surface_idx": {
                                            "type": "integer",
                                            "description": "Volume-local index of the surface the slab belongs to",
                                            "minimum": 0,
                                        },
                                        "thickness": {
                                            "type": "number",
                                            "description": "Radius of the rod",
                                            "minimum": 0,
                                        },
                                        "index_in_coll": {
                                            "type": "integer",
                                            "description": "Volume-local sorting index",
                                            "minimum": 0,
                                        },
                                        "material": {
                                            "type": "object",
                                            "properties": {
                                                "params": {
                                                    "type": "array",
                                                    "description": "Material parameters",
                                                    "minitems": 7,
                                                    "maxItems": 7,
                                                    "items": {"type": "number"},
                                                }
                                            },
                                            "required": ["params"],
                                        },
                                    },
                                    "required": [
                                        "type",
                                        "surface_idx",
                                        "thickness",
                                        "material",
                                    ],
                                },
                            },
                        },
                        "required": ["volume_link", "material_slabs"],
                    },
                }
            },
            "required": ["volumes"],
        },
    },
    "required": ["header", "data"],
}
