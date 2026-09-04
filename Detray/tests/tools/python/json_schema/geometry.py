# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Schema of a geometry file"""

geometry_schema = {
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
                            "enum": ["geometry"],
                        },
                    },
                    "required": ["version", "detector", "tag"],
                },
                "volume_count": {
                    "type": "integer",
                    "description": "Number of volumes in the detector",
                    "minimum": 0,
                },
                "surface_count": {
                    "type": "integer",
                    "description": "Number of surfaces in the detector",
                    "minimum": 0,
                },
            },
            "required": ["common"],
        },
        "data": {
            "type": "object",
            "description": "Geometry data section (volume-by-volume)",
            "properties": {
                "volumes": {
                    "type": "array",
                    "minitems": 1,
                    "items": {
                        "type": "object",
                        "description": "Data of a single volume",
                        "properties": {
                            "name": {
                                "type": "string",
                                "description": "The volume name",
                            },
                            "index": {
                                "type": "integer",
                                "description": "The volume index",
                                "minimum": 0,
                            },
                            "type": {
                                "type": "integer",
                                "description": "The volume type id (e.g. cylinder)",
                                "minimum": 0,
                                "maximum": 4,
                            },
                            "transform": {
                                "type": "object",
                                "description": "Data of the volume placement matrix",
                                "properties": {
                                    "translation": {
                                        "type": "array",
                                        "description": "Translation vector",
                                        "minitems": 3,
                                        "maxItems": 3,
                                        "items": {"type": "number"},
                                    },
                                    "rotation": {
                                        "type": "array",
                                        "description": "Rotation matrix",
                                        "minitems": 9,
                                        "maxItems": 9,
                                        "items": {"type": "number"},
                                    },
                                },
                                "required": ["translation", "rotation"],
                            },
                            "acc_links": {
                                "type": "array",
                                "description": "Typed indices for the acceleration structures",
                                "minitems": 0,
                                "maxItems": 5,
                                "uniqueItems": True,
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "type": {
                                            "type": "integer",
                                            "description": "Type ID of the acceleration structure",
                                            "minimum": 1,
                                            "maximum": 5,
                                        },
                                        "index": {
                                            "type": "integer",
                                            "description": "Global index of the acceleration structure",
                                            "minimum": 0,
                                        },
                                    },
                                    "required": ["type", "index"],
                                },
                            },
                            "surfaces": {
                                "type": "array",
                                "description": "The contained surfaces",
                                "minitems": 1,
                                "uniqueItems": True,
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "identifier": {
                                            "type": "integer",
                                            "description": "detray internal surface hash",
                                            "minimum": 0,
                                        },
                                        "type": {
                                            "type": "integer",
                                            "description": "Surface type (portal: 0, sensitive: 1, passive: 2)",
                                            "minimum": 0,
                                            "maximum": 2,
                                        },
                                        "source": {
                                            "type": "integer",
                                            "description": "bits that link the detray surface to its source",
                                            "minimum": 0,
                                        },
                                        "index_in_coll": {
                                            "type": "integer",
                                            "description": "Volume-local sorting index",
                                            "minimum": 0,
                                        },
                                        "transform": {
                                            "type": "object",
                                            "description": "Data of the surface placement matrix",
                                            "properties": {
                                                "translation": {
                                                    "type": "array",
                                                    "description": "Translation vector",
                                                    "minitems": 3,
                                                    "maxItems": 3,
                                                    "items": {"type": "number"},
                                                },
                                                "rotation": {
                                                    "type": "array",
                                                    "description": "Rotation matrix",
                                                    "minitems": 9,
                                                    "maxItems": 9,
                                                    "items": {"type": "number"},
                                                },
                                            },
                                            "required": ["translation", "rotation"],
                                        },
                                        "mask": {
                                            "type": "object",
                                            "description": "Surface mask data",
                                            "properties": {
                                                "shape": {
                                                    "type": "integer",
                                                    "description": "Shape ID",
                                                    "minimum": 0,
                                                    "maximum": 12,
                                                },
                                                "volume_link": {
                                                    "type": "integer",
                                                    "description": "Navigation link",
                                                    "minimum": 0,
                                                },
                                                "boundaries": {
                                                    "type": "array",
                                                    "description": "Mask boundary values",
                                                    "minitems": 1,
                                                    "maxItems": 7,
                                                    "items": {"type": "number"},
                                                },
                                            },
                                            "required": [
                                                "shape",
                                                "volume_link",
                                                "boundaries",
                                            ],
                                        },
                                        "material": {
                                            "type": "object",
                                            "description": "Surface material link",
                                            "properties": {
                                                "type": {
                                                    "type": "integer",
                                                    "description": "Material type ID",
                                                    "minimum": 0,
                                                    "maximum": 10,
                                                },
                                                "index": {
                                                    "type": "integer",
                                                    "description": "Index of material instance",
                                                    "minimum": 0,
                                                },
                                            },
                                            "required": ["type", "index"],
                                        },
                                    },
                                    "required": ["type", "source", "transform", "mask"],
                                },
                            },
                        },
                        "required": ["name", "index", "type", "transform", "surfaces"],
                    },
                },
                "volume_grid": {
                    "type": "object",
                    "description": "Data for the volume grid (temporary placeholder)",
                    "properties": {
                        "volume_link": {
                            "type": "integer",
                            "description": "Volume link",
                            "minimum": 0,
                        },
                        "acc_link": {
                            "type": "object",
                            "description": "Typed index for the grids",
                            "properties": {
                                "type": {
                                    "type": "integer",
                                    "description": "Type ID of the acceleration structure",
                                },
                                "index": {
                                    "type": "integer",
                                    "description": "Global index of the acceleration structure",
                                    "minimum": 0,
                                },
                            },
                            "required": ["type", "index"],
                        },
                        "axes": {"type": ["null", "array"]},
                        "bins": {"type": ["null", "array"]},
                    },
                    "required": ["volume_link", "acc_link", "axes", "bins"],
                },
            },
            "required": ["volumes"],
        },
    },
    "required": ["header", "data"],
}
