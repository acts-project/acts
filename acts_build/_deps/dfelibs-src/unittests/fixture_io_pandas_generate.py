#!/usr/bin/env python
# SPDX-License-Identifier: MIT
#
# Generate csv/tsv files with numpy/pandas to be read by the `dfe_nameduple`.

from __future__ import print_function

import os.path

import numpy
import pandas

# order is intentionally different from the C++ definition
DTYPE_RECORD = numpy.dtype([
    ('z', 'i8'),
    ('y', 'i4'),
    ('x', 'i2'),
    # should be ? for boolean; pandas writes it out as True/False which does
    # can not be read in by std::iostream into bool again.
    ('d', 'u1'),
    ('c', 'f8'),
    ('b', 'f4'),
    ('a', 'u8'),
])
DTYPE_ADDITIONAL = numpy.dtype([
    ('f0', 'i4'),
    ('f1', 'i4'),
    ('f2', 'i8'),
])

NRECORDS = 32

def make_records(n):
    """
    Build expected records array.
    """
    r = numpy.empty((n,), dtype=DTYPE_RECORD)
    i = numpy.arange(n)
    r['x'] = i
    r['y'] = -2 * i
    r['z'] = 4 * i
    r['a'] = 8 * i
    r['b'] = 0.23126121 # compute product w/ float not double
    r['b'] *= i
    r['c'] = -42.53425 * i
    r['d'] = ((i % 2) != 0)
    return pandas.DataFrame(r)

def make_extra(n):
    """
    Build extra columns not part of the shared record.
    """
    r = numpy.empty((n,), dtype=DTYPE_ADDITIONAL)
    i = numpy.arange(n)
    r['f0'] = i
    r['f1'] = i
    r['f2'] = i
    return pandas.DataFrame(r)

def generate(datadir='.'):
    def make_path(name):
        return os.path.join(datadir, 'namedtuple-{}'.format(name))

    kw = {
        'float_format': '%.18g',
        'index': False,
    }
    records = make_records(NRECORDS)
    # write records with reordered columns
    records.to_csv(make_path('reordered_columns.csv'), sep=',', **kw)
    records.to_csv(make_path('reordered_columns.tsv'), sep='\t', **kw)
    # write records without some columns
    records_minus = records.drop(['x', 'b'], axis=1)
    records_minus.to_csv(make_path('missing_columns.csv'), sep=',', **kw)
    records_minus.to_csv(make_path('missing_columns.tsv'), sep='\t', **kw)
    # write records with extra columns
    records_plus = pandas.concat([records, make_extra(NRECORDS)], axis=1)
    # reorder columns a bit so the additional fields are interleaved
    records_plus = records_plus[['f0', 'x', 'z', 'f1', 'y', 'a', 'c', 'b', 'f2', 'd']]
    records_plus.to_csv(make_path('extra_columns.csv'), sep=',', **kw)
    records_plus.to_csv(make_path('extra_columns.tsv'), sep='\t', **kw)

if __name__ == '__main__':
    import sys
    generate(*sys.argv[1:])
