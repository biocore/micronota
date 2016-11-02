# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import mkstemp

from skbio import DNA, read, write

from micronota.workflow import validate_seq


class Tests(TestCase):
    def setUp(self):
        _, self.ifile = mkstemp()
        _, self.ofile = mkstemp()

    def test_validate_seq(self):
        seqs = [DNA('ATGC', {'id': 'a', 'description': 'A'}),
                DNA('A', {'id': 'b', 'description': 'B'})]
        with open(self.ifile, 'w') as o:
            for seq in seqs:
                write(seq, into=o, format='fasta')
        for l in [0, 2, 4, 5]:
            validate_seq(self.ifile, 'fasta', l, self.ofile)
            obs = list(read(self.ofile, constructor=DNA, format='fasta'))
            exp = [i for i in seqs if len(i) >= l]
            self.assertEqual(exp, obs)

    def test_validate_seq_duplicate_ids(self):
        seqs = [DNA('A', {'id': 'a'}),
                DNA('T', {'id': 'a'})]
        with open(self.ifile, 'w') as o:
            for seq in seqs:
                write(seq, into=o, format='fasta')

        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            validate_seq(self.ifile, 'fasta', 1, self.ofile)


if __name__ == '__main__':
    main()
