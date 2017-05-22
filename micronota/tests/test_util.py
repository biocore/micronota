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
from shutil import rmtree

from skbio import write, read, Sequence
from skbio.metadata import IntervalMetadata

from micronota.util import _filter_sequence_ids, filter_partial_genes


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

    def test_filter_partial_genes(self):
        in_fp = join(self.tmpd, 'in.gff')
        out_fp = join(self.tmpd, 'out.gff')
        imd1 = IntervalMetadata(None)
        imd1.add([(0, 100)], metadata={'partial': '01', 'phase': 0, 'source': 'Prodigal_v2.6.3', 'strand': '.', 'type': '.', 'score': '.'})
        imd2 = IntervalMetadata(None)
        imd2.add([(200, 300)], metadata={'partial': '10', 'phase': 1, 'source': 'Prodigal_v2.6.3', 'strand': '-', 'type': 'CDS', 'score': '1'})
        imd2.add([(2000, 3000)], metadata={'partial': '00', 'phase': 1, 'source': 'Prodigal_v2.6.3', 'strand': '.', 'type': '.', 'score': '.'})

        imd3 = IntervalMetadata(None)
        imd3.add([(2000, 3000)], metadata={'partial': '00', 'phase': 1, 'source': 'Prodigal_v2.6.3', 'strand': '.', 'type': '.', 'score': '.'})

        data = (('seq1', imd1), ('seq2', imd2))
        write(((sid, imd) for sid, imd in data), into=in_fp, format='gff3')
        filter_partial_genes(in_fp, out_fp)
        obs = read(out_fp, format='gff3')
        for i, j in zip(obs, [('seq2', imd3)]):
            self.assertEqual(i, j)

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
