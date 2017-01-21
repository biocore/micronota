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

from pkg_resources import resource_filename
from skbio import DNA, read, write
import yaml

from micronota.workflow import validate_seq, annotate


class Tests(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()
        self.i = join(self.tmpd, 'i.fna')
        self.o = join(self.tmpd, 'o.fna')

    def test_validate_seq(self):
        seqs = [DNA('ATGC', {'id': 'a', 'description': 'A'}),
                DNA('A', {'id': 'b', 'description': 'B'})]
        write((seq for seq in seqs), into=self.i, format='fasta')
        for l in [0, 2, 4, 5]:
            out = '%s_%d' % (self.o, l)
            validate_seq(self.i, 'fasta', l, out)
            obs = list(read(out, constructor=DNA, format='fasta'))
            exp = [i for i in seqs if len(i) >= l]
            self.assertEqual(exp, obs)

    def test_validate_seq_duplicate_ids(self):
        seqs = [DNA('A', {'id': 'a'}),
                DNA('T', {'id': 'a'})]
        write((seq for seq in seqs), into=self.i, format='fasta')

        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            validate_seq(self.i, 'fasta', 0, self.o)

    def test_validate_seq_genbank(self):
        write(DNA('ATGC',
                  {'LOCUS': {'date': '26-APR-1993',
                             'division': 'BCT',
                             'locus_name': 'ECOALKP',
                             'mol_type': 'mRNA',
                             'shape': 'linear',
                             'size': 4,
                             'unit': 'bp'}}),
              into=self.i, format='genbank')
        validate_seq(self.i, 'genbank', 2, self.o)
        exp = '>ECOALKP\nATGC\n'
        with open(self.o) as f:
            obs = f.read()
        self.assertEqual(exp, obs)

    def test_annotate(self):
        config = {'CRISPR': {'minced':
                             {'params': '',
                              'priority': 50,
                              'threads': 1}},
                  'gene': {'prodigal':
                           {'params': '-p meta',
                            'priority': 90,
                            'threads': 1}},
                  'general' : {'metadata': 'foo.sqlite'}}
        config_fp = join(self.tmpd, 'config.yaml')
        with open(config_fp, 'w') as f:
            yaml.dump(config, f, default_flow_style=True)
        write(DNA('ATGC', {'id': 'seq1'}), into=self.i, format='fasta')
        print(self.tmpd)
        annotate(self.i, 'fasta', 1, self.tmpd, 'gff3', 11, 1, True, False, config_fp)
        output = join(self.tmpd, splitext(self.i)[0] + '.valid.fna')
        self.assertTrue(exists(output))
        self.assertTrue(exists(output + '.gff'))

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
