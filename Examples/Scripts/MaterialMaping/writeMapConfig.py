# This file is part of the Acts project.
#
# Copyright (C) 2020 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json
import sys

# Should be run with Python 3 if possible
# Script that parse a Json surfaces map to create an easy to use json config file for the mapping
# Take two arguments in input : The path to the surfaces map and the path of the json config file
# By default the input is : 'surfaces-map.json' and the output is : 'config-map.json'
# The config file can be used to define a binning for all the surfaces in a given volume
# It can also be used to define the binning for volume mapping

def getSurfaceMateral(mat):
    outputmat = {}
    value = {}
    outputmat['Volume'] = mat['Volume']
    if '_Boundary' in mat:
        outputmat['_Boundary'] = 'X'
    if '_Layer' in mat:
        if '__Approach' not in entry:
            if '__Sensitive' not in entry:        
                outputmat['_Layer'] = 'X'
    if '__Approach' in mat:
        outputmat['__Approach'] = mat['__Approach']
    if '__Sensitive' in mat:
        outputmat['_Layer'] = mat['_Layer']
        outputmat['__Sensitive'] = 'X'
    value['binUtility'] = mat['value']['binUtility']
    value['mapMaterial'] = False
    value['stype'] = mat['value']['stype']
    outputmat['value'] = value
    return outputmat

if sys.version_info[0] < 3:
    print('Using Python 2')
    print('To obtain the proper ordering in the Json files Python 3 is recomanded')

if len(sys.argv) < 2 :
    inFileName = 'geometry-maps.json'
else :
    inFileName = sys.argv[1]

    
with open(inFileName,'r') as json_file:
    config = {}
    data = json.load(json_file)
    lastVol = -1
    for entry in data['2.Surfaces']['entries']:
        if lastVol != entry['Volume']:
            if lastVol != -1:
                config[lastVol] = vconfig
            vconfig = []
            lastVol = entry['Volume']
            typeLayer = []
            createdApproach1 = False
            createdApproach2 = False
            typeBoundary = []
            typeSensitive = []
            listLayer = []

        if 'stype' not in entry['value']:
            entry['value']['stype'] = ''

        if '_Layer' in entry:  
            if '__Approach' not in entry:
                if '__Sensitive' not in entry:
                    if entry['value']['stype'] not in typeLayer:
                        typeLayer.append(entry['value']['stype'])
                        surface = getSurfaceMateral(entry)
                        vconfig.append(surface)
                        continue

        if '_Boundary' in entry:    
            if '_Layer' not in entry:
                if entry['value']['stype'] not in typeBoundary:
                    typeBoundary.append(entry['value']['stype'])
                    surface = getSurfaceMateral(entry)
                    vconfig.append(surface)
                    continue         

        if '__Approach' in entry:
            if '__Sensitive' not in entry:
                if entry['__Approach'] == 1 and createdApproach1 == False:
                    createdApproach1 = True
                    surface = getSurfaceMateral(entry)
                    vconfig.append(surface)
                    continue
                if entry['__Approach'] == 2 and createdApproach2 == False:
                    createdApproach2 = True
                    surface = getSurfaceMateral(entry)
                    vconfig.append(surface)
                    continue

        if '__Sensitive' in entry:  
            if '__Approach' not in entry:
                if entry['value']['stype'] not in typeSensitive:
                    if entry['_Layer'] not in listLayer:
                        listLayer.append(entry['_Layer'])
                        typeSensitive.append(entry['value']['stype'])
                        surface = getSurfaceMateral(entry)
                        vconfig.append(surface)
                        continue

    if lastVol != -1:
        config[lastVol] = vconfig
    config['2.Surfaces'] = vconfig

if len(sys.argv) < 3 :
    outFileName = 'config-map.json'
else :
    outFileName = sys.argv[2]
    
with open(outFileName, 'w') as outfile:
    json.dump(config, outfile, indent=4)
