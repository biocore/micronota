# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkdtemp
from os.path import join, splitext, exists
from shutil import rmtree
from io import StringIO

import yaml
from skbio import DNA, read, write
from skbio.metadata import IntervalMetadata
from skbio.util import get_data_path

from micronota.workflow import annotate, summarize, create_faa


class Tests(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()
        self.i = join(self.tmpd, 'i.fna')
        self.o = join(self.tmpd, 'o.fna')

    def test_annotate(self):
        config = {'structural_annotation': {'minced':
                                            {'params': '',
                                             'priority': 50,
                                             'output': 'minced',
                                             'threads': 1},
                                            'prodigal':
                                            {'params': '-p meta -f gff',
                                             'priority': 90,
                                             'output': 'prodigal',
                                             'threads': 1}},
                  'protein': {},
                  'bacteria': {},
                  'general': {'metadata': 'foo.sqlite'}}
        config_fp = join(self.tmpd, 'config.yaml')
        with open(config_fp, 'w') as f:
            yaml.dump(config, f, default_flow_style=True)
        write(DNA('ATGC', {'id': 'seq1'}), into=self.i, format='fasta')
        annotate(self.i, 'fasta', 1, self.tmpd, 'gff3', 11, 'bacteria', 'metagenome', (), 1,
                 True, False, True, config_fp)
        output = join(self.tmpd, splitext(self.i)[0])
        self.assertTrue(exists(output + '.fna'))
        self.assertTrue(exists(output + '.gff3'))

    def test_summarize(self):
        gff = get_data_path('summarize.gff')
        seqs = [DNA('A' * 5000000, metadata={'id': 'gi|556503834|ref|NC_000913.3|'}),
                DNA('AG' * 2500000, metadata={'id': 'gi|556503834|ref|NC_000913.2|'})]
        for (seq_id, imd), seq in zip(read(gff, format='gff3'), seqs):
            seq.interval_metadata = imd

        with StringIO() as obs, open(get_data_path('summarize.txt')) as exp:
            summarize(seqs, obs)
            self.assertEqual(obs.getvalue(), exp.read())

    def test_create_faa(self):
        imd = IntervalMetadata(None)
        imd.add([(0, 120)], metadata={'type': 'CDS', 'product': 'Homoserine kinase',
                                      'ID': '1_1'})
        seq = DNA('ATGGTTAAAGTTTATGCCCCGGCTTCCAGTGCCAATATGAGCGTCGGGTTTGATGTGCTC'
                  'GGGGCGGCGGTGACACCCGTTGATGGTGCATTGCTCGGAGATGTAGTCACGGTTGAGGCG'
                  'GCAGAGACATTCAGTCTCAACAACCTCGGACGCTTTGCCGATAAGCTGCCGTCAGAACCA'
                  'CGGGAAAATATCGTTTATCAGTGCTGGGAGCGTTTTTGCCTGGAGCTGGGCAAGCAAATT'
                  'CCAGTGGCGATGACTCTGGAAAAGAATATGCCGATCGGCTCGGGCTTAGGCTCCAGCGCC'
                  'TGTTCGGTGGTCGCGGCGCTGATGGCGATGAATGAACACTGCGGCAAGCCACTTAATGAC'
                  'ACCCGTTTGCTGGCTTTGATGGGCGAGCTGGAAGGACGTATCTCCGGCAGCATTCATTAC'
                  'GACAACGTGGCACCGTGTTTTCTTGGTGGTATGCAGTTGATGATCGAAGAAAACGACATC'
                  'ATCAGCCAGCAAGTGCCAGGGTTTGATGAGTGGCTGTGGGTGCTGGCGTATCCGGGAATT'
                  'AAAGTCTCGACGGCAGAAGCCCGGGCTATTTTACCGGCGCAGTATCGCCGCCAGGATTGC'
                  'ATTGCGCACGGGCGACATCTGGCTGGCTTCATTCACGCCTGCTATTCCCGTCAGCCTGAG'
                  'CTTGCCGCGAAGCTGATGAAAGATGTTATCGCTGAACCCTACCGTGAACGGTTACTGCCT'
                  'GGCTTCCGGCAGGCGCGGCAGGCGGTCGCGGAAATCGGCGCGGTAGCGAGCGGTATCTCC'
                  'GGCTCCGGCCCGACCTTGTTCGCTCTATGTGACAAGCCGGATACCGCCCAGCGCGTTGCC'
                  'GACTGGTTGGGTAAGAACTACCTGCAAAATCAGGAAGGTTTTGTTCATATTTGCCGGCTG'
                  'GATACGGCGGGCGCACGAGTACTGGAAAACTAA',
                  interval_metadata=imd)
        create_faa([seq], self.o)
        exp = '>1_1 Homoserine kinase\nMVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEA\n'

        with open(self.o) as out:
            obs = out.read()
            self.assertEqual(exp, obs)

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
