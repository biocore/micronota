# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, mock
from os.path import join, dirname, relpath
from tempfile import mkdtemp
from shutil import rmtree

from skbio import read, write

from micronota.util import _create_config
from micronota.workflow import annotate


class TestAnnotate(TestCase):
    def setUp(self):
        test_dir = join('micronota', 'db', 'tests', 'data', 'uniref')
        cwd = dirname(__file__)
        test_dir = '/Users/zech/Dropbox/git/micronota/micronota/db/tests/data/uniref'
        files = [
            'uniref100_Swiss-Prot_Archaea.fna',
            'uniref100_Swiss-Prot_Bacteria.fna',
            'uniref100_Swiss-Prot_Eukaryota.fna',
            'uniref100_Swiss-Prot_Viruses.fna',
            'uniref100_TrEMBL_Archaea.fna',
            'uniref100_TrEMBL_Bacteria.fna',
            'uniref100_TrEMBL_Eukaryota.fna',
            'uniref100_TrEMBL_Viruses.fna']
        files = [join(test_dir, f) for f in files]
        self.tmp = mkdtemp()
        self.test1 = join(self.tmp, 'test1.fna')
        with open(self.test1, 'w') as f:
            for seq in read(files[1], format='fasta'):
                write(seq, format='fasta', into=f)

        self.obs_tmp = mkdtemp()
        # mock up the db path
        self.patcher = mock.patch('micronota.workflow._DB_PATH', test_dir)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()
        rmtree(self.tmp)
        rmtree(self.obs_tmp)

    def test_annotate(self):
        config = _create_config(None)
        annotate(self.test1, 'fasta', self.obs_tmp, 'genbank', 'archaea', 1, config)


if __name__ == '__main__':
    main()
