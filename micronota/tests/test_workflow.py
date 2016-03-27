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
from filecmp import cmp

from skbio import read, write
from skbio.util import get_data_path
from skbio.metadata import Feature

from micronota.bfillings.diamond import DiamondCache
from micronota.workflow import annotate, annotate_all_cds
from micronota.config import Configuration


class TestAnnotate(TestCase):
    def setUp(self):
        self.test_dir = abspath(
            join('micronota', 'db', 'tests', 'data', 'uniref', 'uniref100'))
        files = [
            'Swiss-Prot_Archaea.fna',
            'Swiss-Prot_Bacteria.fna',
            'Swiss-Prot_Eukaryota.fna',
            'Swiss-Prot_Viruses.fna',
            'TrEMBL_Archaea.fna',
            'TrEMBL_Bacteria.fna',
            'TrEMBL_Eukaryota.fna',
            'TrEMBL_Viruses.fna']
        files = [join(self.test_dir, f) for f in files]
        self.tmp = mkdtemp()
        self.test1 = join(self.tmp, 'test1.fna')
        self.test1_exp = 'test1.genbank'
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
                 1, 'archaea', True, config, cache=False)
        self.assertTrue(cmp(
            get_data_path(self.test1_exp),
            join(self.obs_tmp, self.test1_exp),
            shallow=False))

    def test_annotate_cache(self):
        config = Configuration()
        config.db_dir = self.test_dir
        annotate(self.test1, 'fasta', self.obs_tmp, 'genbank',
                 1, 'archaea', True, config, cache=True)
        self.assertTrue(cmp(
            get_data_path(self.test1_exp),
            join(self.obs_tmp, self.test1_exp),
            shallow=False))


class TestAnnotateCDS(TestCase):
    def setUp(self):
        seq = ('MNSFRKTCAGALALIFGATSIVPTVAAPMNMDRPAINQNVIQARAHYR'
               'PQNYNRGHRPGYWHGHRGYRHYRHGYRRHNDGWWYPLAAFGAGAIIGG'
               'AISQPRPVYRAPAGSPHVQWCYSRYKSYRASDNTFQPYNGPRKQCRSP'
               'YSR')

        note = "start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.678"
        self.test_dir = abspath(
            join('micronota', 'db', 'tests', 'data', 'uniref', 'uniref100'))
        self.im = {
            Feature(note=note,
                    id="1_1",
                    translation=seq,
                    location='<1..>441',
                    left_partial_=True,
                    rc_=False,
                    type_='CDS',
                    right_partial_=True): [(0, 441)]}
        self.obs_tmp = mkdtemp()

    def test_annotate_all_cds_no_cache(self):

        config = Configuration()
        config.db_dir = self.test_dir
        im, cache = annotate_all_cds(self.im, out_dir=self.obs_tmp,
                                     kingdom='archaea', config=config,
                                     cpus=1, cache=None)
        self.assertTrue(cache is None)

    def test_annotate_all_cds_cache_enabled(self):

        config = Configuration()
        config.db_dir = self.test_dir
        im, cache = annotate_all_cds(self.im, out_dir=self.obs_tmp,
                                     kingdom='archaea', config=config,
                                     cpus=1, cache=DiamondCache())
        self.assertTrue(cache is not None)
        self.assertTrue(len(cache.seqs) > 0)


if __name__ == '__main__':
    main()
