#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from tempfile import mkstemp
from os import getcwd, remove, close
from os.path import join
from functools import partial
from unittest import TestCase, main

from skbio.util import get_data_path
from burrito.util import ApplicationError

from micronota.bfillings.infernal import (
    CMScan, cmscan_fasta,
    cmpress_cm)


class InfernalTests(TestCase):
    def setUp(self):
        self.temp_fd, self.temp_fp = mkstemp()
        self.get_infernal_path = partial(
            get_data_path, subfolder=join('data', 'infernal'))

        self.positive_fps = list(map(self.get_infernal_path, [
            # modified from NC_018498.gbk
            'NC_018498.fna',
            ]))
        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only']))
        self.cm_fp = self.get_infernal_path('RF00522.cm')

    def tearDown(self):
        # remove the tempdir and contents
        close(self.temp_fd)
        remove(self.temp_fp)


class CMScanTests(InfernalTests):
    def test_base_command(self):
        c = CMScan()
        self.assertEqual(
            c.BaseCommand,
            'cd "%s/"; %s' % (getcwd(), c._command))
        params = {'--rfam': None, '--cpu': 2}
        for i in params:
            if params[i] is None:
                c.Parameters[i].on()
                cmd = 'cd "{d}/"; {cmd} {option}'.format(
                    d=getcwd(), cmd=c._command, option=i)
            else:
                c.Parameters[i].on(params[i])
                cmd = 'cd "{d}/"; {cmd} {option} {value}'.format(
                    d=getcwd(), cmd=c._command,
                    option=i, value=params[i])

            self.assertEqual(c.BaseCommand, cmd)
            c.Parameters[i].off()

    def test_cmscan_fasta_wrong_input(self):
        for fp in self.negative_fps:
            with self.assertRaisesRegex(
                    ApplicationError,
                    r'Error: Sequence file .* is empty or misformatted'):
                cmscan_fasta(self.cm_fp, fp, 'foo')

    def test_cmscan_fasta(self):
        params = {'--rfam': None, '--noali': None}
        for f in self.positive_fps:
            res = cmscan_fasta(self.cm_fp, f, self.temp_fp, 0.1, 1, params)
            res['StdOut'].close()
            res['StdErr'].close()
            obs = res['--tblout']
            out_fp = '.'.join([f, 'tblout'])
            with open(out_fp) as exp:
                # skip comment lines as some contain running time info
                self.assertListEqual(
                    [i for i in exp.readlines() if not i.startswith('#')],
                    [j for j in obs.readlines() if not j.startswith('#')])
            obs.close()

    def test_cmscan_sequence(self):
        pass


class CMPressTests(InfernalTests):
    def test_cmpress_cm_exist(self):
        with self.assertRaisesRegex(
                ApplicationError,
                r'Error: Looks like .* is already pressed'):
            cmpress_cm(self.cm_fp)

    def test_compress_cm(self):
        # .i1i file is different from run to run. skip it.
        suffices = ('i1f', 'i1m', 'i1p')
        exp = []
        for i in suffices:
            with open('.'.join([self.cm_fp, i]), 'rb') as f:
                exp.append(f.read())

        res = cmpress_cm(self.cm_fp, True)
        res['StdOut'].close()
        res['StdErr'].close()

        for i, e in zip(suffices, exp):
            with open('.'.join([self.cm_fp, i]), 'rb') as f:
                self.assertEqual(f.read(), e)


if __name__ == '__main__':
    main()
