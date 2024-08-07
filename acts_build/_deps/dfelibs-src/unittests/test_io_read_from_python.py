#!/usr/bin/env python
# SPDX-License-Identifier: MIT
#
# Read in files written by `dfe_namedtuple.hpp` with numpy and check content

from __future__ import print_function

import os.path

import numpy

DTYPE = numpy.dtype([
    ('x', 'i2'),
    ('y', 'i4'),
    ('z', 'i8'),
    ('a', 'u8'),
    ('b', 'f4'),
    ('c', 'f8'),
    ('d', '?'),
])
NRECORDS = 1024

def make_records(n):
    """
    Build expected records array.
    """
    r = numpy.empty((n,), dtype=DTYPE)
    i = numpy.arange(n)
    r['x'] = i
    r['y'] = -2 * i
    r['z'] = 4 * i
    r['a'] = 8 * i
    r['b'] = 0.23126121 # compute product w/ float not double
    r['b'] *= i
    r['c'] = -42.53425 * i
    r['d'] = ((i % 2) != 0)
    return r

def _check(name, records):
    """
    Check that the given records array matches the expected.
    """
    if records.dtype != DTYPE:
        raise TypeError('inconsistent dtypes', name)
    compare = make_records(NRECORDS)
    for column in DTYPE.names:
        matches = (records[column] == compare[column])
        if not matches.all():
            raise ValueError('inconsistent column value', name, column)
        print('{} column {} matches'.format(name, column))

def test(datadir='.'):
    _check('csv', numpy.loadtxt(os.path.join(datadir, 'test.csv'),
        dtype=DTYPE, delimiter=',', skiprows=1))
    _check('tsv', numpy.loadtxt(os.path.join(datadir, 'test.tsv'),
        dtype=DTYPE, delimiter='\t', skiprows=1))
    _check('npy', numpy.load(os.path.join(datadir, 'test.npy')))

if __name__ == '__main__':
    import sys
    test(*sys.argv[1:])
