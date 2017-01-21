# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join
from tempfile import mkdtemp

from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata
from skbio import write, read, Sequence

from micronota.util import _filter_sequence_ids


class Tests(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()

    def test_filter_sequence_ids(self):
        seqs = [Sequence('A', {'id': 'seq1', 'description': ''}),
                Sequence('T', {'id': 'seq2', 'description': ''})]

        ifile = join(self.tmpd, 'in.fna')
        write((i for i in seqs), into=ifile, format='fasta')
        ofile = join(self.tmpd, 'out.fna')

        idss = [('foo'), {'seq1'}, ('seq2'), {'seq1', 'seq2'}]
        exps = [seqs, seqs[1:], seqs[:-1], []]

        for ids, exp in zip(idss, exps):
            _filter_sequence_ids(ifile, ofile, ids)
            obs = list(read(ofile, constructor=Sequence, format='fasta'))
            self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()

