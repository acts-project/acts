# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Schema of a surface grid file"""

surface_grid_schema = {
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
                            "enum": ["surface_grids"],
                        },
                    },
                    "required": ["version", "detector", "tag"],
                },
                "grid_count": {
                    "type": "integer",
                    "description": "Number of grids in the detector",
                    "minimum": 0,
                },
            },
            "required": ["common"],
        },
        "data": {
            "type": "object",
            "description": "Grid data section (volume-by-volume)",
            "properties": {
                "grids": {
                    "type": "array",
                    "minitems": 1,
                    "items": {
                        "type": "object",
                        "description": "Data of all grids in a single volume",
                        "properties": {
                            "volume_link": {
                                "type": "integer",
                                "description": "Volume link",
                                "minimum": 0,
                            },
                            "grid_data": {
                                "type": "array",
                                "minitems": 1,
                                "items": {
                                    "type": "object",
                                    "description": "Data of a single grid",
                                    "properties": {
                                        "owner_link": {
                                            "type": "integer",
                                            "description": "Link to the owner of the grid",
                                            "minimum": 0,
                                        },
                                        "grid_link": {
                                            "type": "object",
                                            "description": "Typed index for the grid",
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
                                        "axes": {
                                            "type": "array",
                                            "descriptrion": "The axes of the grid",
                                            "minitems": 1,
                                            "uniqueItems": True,
                                            "items": {
                                                "type": "object",
                                                "properties": {
                                                    "label": {
                                                        "type": "integer",
                                                        "description": "Label ID of the axis",
                                                        "minimum": 0,
                                                        "maximum": 2,
                                                    },
                                                    "bounds": {
                                                        "type": "integer",
                                                        "description": "Type ID of the bounds",
                                                        "minimum": 0,
                                                        "maximum": 2,
                                                    },
                                                    "binning": {
                                                        "type": "integer",
                                                        "description": "Type ID of the binning",
                                                        "minimum": 0,
                                                        "maximum": 1,
                                                    },
                                                    "bins": {
                                                        "type": "integer",
                                                        "description": "Number of bins",
                                                        "minimum": 1,
                                                    },
                                                    "edges": {
                                                        "type": "array",
                                                        "description": "Bin edges (min, max for regular binning)",
                                                        "minitems": 2,
                                                        "uniqueItems": True,
                                                        "items": {"type": "number"},
                                                    },
                                                },
                                                "required": [
                                                    "label",
                                                    "bounds",
                                                    "binning",
                                                    "bins",
                                                    "edges",
                                                ],
                                            },
                                        },
                                        "bins": {
                                            "type": "array",
                                            "description": "The bin content",
                                            "minitems": 1,
                                            "uniqueItems": True,
                                            "items": {
                                                "type": "object",
                                                "description": "A single bin",
                                                "properties": {
                                                    "loc_index": {
                                                        "type": "array",
                                                        "description": "Local bin indices",
                                                        "minitems": 1,
                                                        "maxitems": 3,
                                                        "items": {
                                                            "type": "integer",
                                                            "minimum": 0,
                                                        },
                                                    },
                                                    "content": {
                                                        "type": "array",
                                                        "description": "Global surface indices",
                                                        "items": {
                                                            "type": "integer",
                                                            "minimum": 0,
                                                        },
                                                    },
                                                },
                                                "required": ["loc_index", "content"],
                                            },
                                        },
                                    },
                                    "required": ["owner_link", "axes", "bins"],
                                },
                            },
                        },
                        "required": ["volume_link", "grid_data"],
                    },
                }
            },
            "required": ["grids"],
        },
    },
    "required": ["header", "data"],
}
