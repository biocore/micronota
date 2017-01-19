r'''
Utility functionality
=====================

.. currentmodule:: micronota.util

This module (:mod:`micronota.util`) provides various utility functionality,
'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import shutil
from os import remove
from os.path import exists, isdir, join, abspath, dirname, basename, splitext
from urllib.request import urlopen
from unittest import TestCase
from sqlite3 import connect
from inspect import stack

from skbio import read, write, Sequence


def _filter_sequence_ids(in_fp, out_fp, ids):
    '''Filter away the seq with specified IDs.'''
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format='fasta', constructor=Sequence):
            seq_id = seq.metadata['id']
            if seq_id not in ids:
                write(seq, format='fasta', into=out)


class _DBTest(TestCase):
    def _test_eq_db(self, db1, db2):
        '''Test if two database files have the same contents.'''
        with connect(db1) as o, connect(db2) as e:
            co = o.cursor()
            ce = e.cursor()
            # compare each table
            for table, in e.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"):
                co.execute('SELECT * FROM %s' % table)
                ce.execute('SELECT * FROM %s' % table)
                self.assertEqual(co.fetchall(), ce.fetchall())
