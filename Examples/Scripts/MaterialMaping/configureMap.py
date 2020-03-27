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
# Script that use the json config file to configure the Json surfaces map for the material mapping
# Take two arguments in input : The path to the surfaces map and the path of the json config file
# By default the input is : 'surfaces-map.json' and the output is : 'config-map.json'
# The config file can be used to define a binning for all the surface of a specific type and boundary type in a given volume
# If the binning is kept to (1,1) the surface will not be mapped on.

if sys.version_info[0] < 3:
    print('Using Python 2')
    print('To obtain the proper ordering in the Json files Python 3 is recomanded')

if len(sys.argv) < 2 :
    inFileName = 'surfaces-map.json'
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

        for kvol in data['volumes']:
            Name = data['volumes'][kvol]['Name']

            if 'boundaries' in data['volumes'][kvol] :
                for kbound in data['volumes'][kvol]['boundaries'] :
                    dbound = data['volumes'][kvol]['boundaries'][kbound]
                    data['volumes'][kvol]['boundaries'][kbound]['bin0'] = config[Name]['boundaries'][dbound['stype']][0]
                    data['volumes'][kvol]['boundaries'][kbound]['bin1'] = config[Name]['boundaries'][dbound['stype']][1]
                    if (config[Name]['boundaries'][dbound['stype']][0][2] != 1) or (config[Name]['boundaries'][dbound['stype']][1][2] != 1) :
                        data['volumes'][kvol]['boundaries'][kbound]['matSurface'] = True

            if 'layers' in data['volumes'][kvol] :
                for klay in data['volumes'][kvol]['layers'] :

                    if 'representing' in data['volumes'][kvol]['layers'][klay] :
                        drep = data['volumes'][kvol]['layers'][klay]['representing']
                        data['volumes'][kvol]['layers'][klay]['representing']['bin0'] = config[Name]['representing'][drep['stype']][0]
                        data['volumes'][kvol]['layers'][klay]['representing']['bin1'] = config[Name]['representing'][drep['stype']][1]
                        if (config[Name]['representing'][drep['stype']][0][2] != 1) or (config[Name]['representing'][drep['stype']][1][2] != 1) :
                             data['volumes'][kvol]['layers'][klay]['representing']['matSurface'] = True

                    if 'approach' in data['volumes'][kvol]['layers'][klay] :
                        for kapp  in data['volumes'][kvol]['layers'][klay]['approach'] :
                            dapp = data['volumes'][kvol]['layers'][klay]['approach'][kapp]
                            data['volumes'][kvol]['layers'][klay]['approach'][kapp]['bin0'] = config[Name]['approach'][dapp['stype']][kapp][0]
                            data['volumes'][kvol]['layers'][klay]['approach'][kapp]['bin1'] = config[Name]['approach'][dapp['stype']][kapp][1]
                            if (config[Name]['approach'][dapp['stype']][kapp][0][2] != 1) or (config[Name]['approach'][dapp['stype']][kapp][1][2] != 1) :
                                data['volumes'][kvol]['layers'][klay]['approach'][kapp]['matSurface'] = True


                    if 'sensitive' in data['volumes'][kvol]['layers'][klay] :
                        for ksen  in data['volumes'][kvol]['layers'][klay]['sensitive'] :
                            dsen = data['volumes'][kvol]['layers'][klay]['sensitive'][ksen]
                            data['volumes'][kvol]['layers'][klay]['sensitive'][ksen]['bin0'] = config[Name]['sensitive'][dsen['stype']][0]
                            data['volumes'][kvol]['layers'][klay]['sensitive'][ksen]['bin1'] = config[Name]['sensitive'][dsen['stype']][1]
                            if (config[Name]['sensitive'][dsen['stype']][0][2] != 1) or (config[Name]['sensitive'][dsen['stype']][1][2] != 1) :
                                data['volumes'][kvol]['layers'][klay]['sensitive'][ksen] = True

    json_file.seek(0)
    json.dump(data, json_file, indent=4)
    json_file.truncate()
