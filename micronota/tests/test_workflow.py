# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join, abspath
from tempfile import mkdtemp
from shutil import rmtree

from skbio import read, write

from micronota.workflow import annotate
from micronota.config import Configuration


class TestAnnotate(TestCase):
    def setUp(self):
        self.test_dir = abspath(
            join('micronota', 'db', 'tests', 'data', 'uniref'))
        files = [
            'uniref100_Swiss-Prot_Archaea.fna',
            'uniref100_Swiss-Prot_Bacteria.fna',
            'uniref100_Swiss-Prot_Eukaryota.fna',
            'uniref100_Swiss-Prot_Viruses.fna',
            'uniref100_TrEMBL_Archaea.fna',
            'uniref100_TrEMBL_Bacteria.fna',
            'uniref100_TrEMBL_Eukaryota.fna',
            'uniref100_TrEMBL_Viruses.fna']
        files = [join(self.test_dir, f) for f in files]
        self.tmp = mkdtemp()
        self.test1 = join(self.tmp, 'test1.fna')
        with open(self.test1, 'w') as f:
            for seq in read(files[1], format='fasta'):
                write(seq, format='fasta', into=f)

        self.obs_tmp = mkdtemp()

    def tearDown(self):
        rmtree(self.tmp)
        rmtree(self.obs_tmp)

    def test_annotate(self):
        config = Configuration()
        config.db_dir = self.test_dir
        annotate(self.test1, 'fasta', self.obs_tmp, 'genbank',
                 1, 'archaea', True, config)


if __name__ == '__main__':
    main()
