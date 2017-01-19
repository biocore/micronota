# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkdtemp
from os.path import join
from shutil import rmtree

from skbio import DNA, read, write

from micronota.workflow import validate_seq


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

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
