# This file is part of the Acts project.
#
# Copyright (C) 2020-2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json
import sys

# Should be run with Python 3 if possible
# Script that use the json config file to configure the Json surfaces map for the material mapping
# Take two arguments in input : The path to the surfaces map and the path of the json config file
# By default the input is : 'surfaces-map.json' and the output is : 'config-map.json'
# The config file can be used to define a binning for all the surfaces in a given volume
# It can also be used to define the binning for volume mapping

if sys.version_info[0] < 3:
    print('Using Python 2')
    print('To obtain the proper ordering in the Json files Python 3 is recomanded')

if len(sys.argv) < 2 :
    inFileName = 'geometry-maps.json'
    confFileName = 'config-map.json'
    
if len(sys.argv) < 3 :
    confFileName = 'config-map.json'
    
else :
    inFileName = sys.argv[1]
    confFileName = sys.argv[2]

    
with open(inFileName,'r+') as json_file:
    with open(confFileName,'r') as config_file:

        config = json.load(config_file)
        data = json.load(json_file)

        for entry in data['2.Surfaces']['entries']:

            if 'stype' not in entry['value']:
                entry['value']['stype'] = ''

            if '_Layer' in entry:  
                if '__Approach' not in entry:
                    if '__Sensitive' not in entry:
                        for conf in config[str(entry['Volume'])]:
                            if '_Layer' in conf and conf['_Layer'] == 'X' and conf['value']['stype'] == entry['value']['stype']:
                                entry['value']['mapMaterial'] = conf['value']['mapMaterial']
                                ibin = 0
                                for bin in entry['value']['binUtility']['binningdata']:                                  
                                    bin['bins'] = conf['value']['binUtility']['binningdata'][ibin]['bins']
                                    ibin = ibin+1
                                continue
                        continue

            if '_Boundary' in entry:    
                if '_Layer' not in entry:
                    for conf in config[str(entry['Volume'])]:
                        if '_Boundary' in conf and conf['_Boundary'] == 'X' and conf['value']['stype'] == entry['value']['stype']:
                            entry['value']['mapMaterial'] = conf['value']['mapMaterial']
                            ibin = 0
                            for bin in entry['value']['binUtility']['binningdata']:
                                bin['bins'] = conf['value']['binUtility']['binningdata'][ibin]['bins']
                                ibin = ibin+1
                            continue
                    continue
                 
            if '__Approach' in entry:
                if '__Sensitive' not in entry:
                    for conf in config[str(entry['Volume'])]:
                        if '__Approach' in conf and conf['__Approach'] == entry['__Approach'] and conf['value']['stype'] == entry['value']['stype']:
                            entry['value']['mapMaterial'] = conf['value']['mapMaterial']
                            ibin = 0
                            for bin in entry['value']['binUtility']['binningdata']:
                                bin['bins'] = conf['value']['binUtility']['binningdata'][ibin]['bins']
                                ibin = ibin+1
                            continue
                    continue
                 
            if '__Sensitive' in entry:  
                if '__Approach' not in entry:
                    for conf in config[str(entry['Volume'])]:
                        if '__Sensitive' in conf and conf['__Sensitive'] == 'X' and conf['_Layer'] == entry['_Layer'] and conf['value']['stype'] == entry['value']['stype']:
                            entry['value']['mapMaterial'] = conf['value']['mapMaterial']
                            ibin = 0
                            for bin in entry['value']['binUtility']['binningdata']:
                                bin['bins'] = conf['value']['binUtility']['binningdata'][ibin]['bins']
                                ibin = ibin+1
                            continue
                    continue  
                     
    json_file.seek(0) 
    json.dump(data, json_file, indent=4)
    json_file.truncate()
    