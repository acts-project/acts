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
# The config file can be used to define a binning for all the surface of a specific type and boundary type in a given volume
# If the binning is kept to (1,1) the surface will not be mapped on.

if sys.version_info[0] < 3:
    print('Using Python 2')
    print('To obtain the proper ordering in the Json files Python 3 is recomanded')

if len(sys.argv) < 2 :
    inFileName = 'surfaces-map.json'
else :
    inFileName = sys.argv[1]

with open(inFileName,'r') as json_file:
    config = {}
    data = json.load(json_file)

    for kvol in data['volumes']:
        vconfig = {}
        bconfig = {}
        rconfig = {}
        aconfig = {}
        abin    = {}
        sconfig = {}

        if 'boundaries' in data['volumes'][kvol] :
            for kbound in data['volumes'][kvol]['boundaries'] :
                dbound = data['volumes'][kvol]['boundaries'][kbound]
                if not dbound['stype'] in bconfig :
                    bconfig[dbound['stype']] = [dbound['bin0'], dbound['bin1']]
            vconfig['boundaries']=bconfig

        if 'layers' in data['volumes'][kvol] :
            for klay in data['volumes'][kvol]['layers'] :

                if 'representing' in data['volumes'][kvol]['layers'][klay] :
                    drep = data['volumes'][kvol]['layers'][klay]['representing']
                    if not drep['stype'] in rconfig :
                        rconfig[drep['stype']] = [drep['bin0'], drep['bin1']]
                    vconfig['representing'] = rconfig

                if 'approach' in data['volumes'][kvol]['layers'][klay] :
                    for kapp  in data['volumes'][kvol]['layers'][klay]['approach'] :
                        dapp = data['volumes'][kvol]['layers'][klay]['approach'][kapp]
                        abin[kapp] = [dapp['bin0'], dapp['bin1']]
                    aconfig[dapp['stype']] = abin
                    vconfig['approach'] = aconfig

                if 'sensitive' in data['volumes'][kvol]['layers'][klay] :
                    for ksen  in data['volumes'][kvol]['layers'][klay]['sensitive'] :
                        dsen = data['volumes'][kvol]['layers'][klay]['sensitive'][ksen]
                        if not dsen['stype'] in sconfig :
                            sconfig[dsen['stype']] = [dsen['bin0'], dsen['bin1']]
                    vconfig['sensitive'] = sconfig
        config[data['volumes'][kvol]['Name']] = vconfig


if len(sys.argv) < 3 :
    outFileName = 'config-map.json'
else :
    outFileName = sys.argv[2]

with open(outFileName, 'w') as outfile:
    json.dump(config, outfile, indent=4)
