# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from tempfile import NamedTemporaryFile

from skbio import DNA
from skbio.io import read, write
from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata

from micronota.quality import (
    compute_seq_score, compute_rrna_score,
    compute_trna_score, compute_gene_score)


class ScoreTests(TestCase):
    def test_compute_seq_score(self):
        seqs = [DNA('A' * 10000)] * 2
        obs = compute_seq_score(seqs, contigs=False)
        self.assertEqual(obs, 1)

        obs = compute_seq_score(seqs, contigs=True)
        self.assertEqual(obs, 2/3)

    def test_compute_seq_score_Ns(self):
        seqs = [DNA('A' * 9990 + 'N' * 10),
                DNA('A' * 9000 + 'N' * 1000)]
        obs = compute_seq_score(seqs)
        exp = 2 / (2 + 0 + 1 + 2)
        self.assertEqual(obs, exp)

    def test_compute_seq_score_bad(self):
        seqs = [DNA('A' * 5000 + 'Y' * 5000),
                DNA('A' * 4000 + 'Y' * 5000 + 'N' * 1000)]
        obs = compute_seq_score(seqs)
        exp = 1 / (1 + 1 + 1 + 1)
        self.assertEqual(obs, exp)

    def test_compute_rrna_score(self):
        imd1 = IntervalMetadata(None)
        imd1.add([(0, 100)], metadata={'type': 'rRNA', 'product': '5s_rRNA'})
        imd1.add([(0, 725)], metadata={'type': 'rRNA', 'product': '16s_rRNA'})
        imd1.add([(0, 12)], metadata={'type': 'tRNA'})
        imd2 = IntervalMetadata(None)
        imd2.add([(0, 10)], metadata={'type': 'rRNA', 'product': '5s_rRNA'})
        imd2.add([(0, 1450)], metadata={'type': 'rRNA', 'product': '16s_rRNA'})
        imd3 = IntervalMetadata(None)
        imd3.add([(0, 2900)], metadata={'type': 'rRNA', 'product': '23s_rRNA'})

        obs = compute_rrna_score([imd1, imd2, imd3])
        exp = 0.1 + sum([0.3] * 3)
        self.assertEqual(obs, exp)

        obs = compute_rrna_score([imd1, imd2])
        exp = 0.1 + sum([0.3] * 2)
        self.assertEqual(obs, exp)

        obs = compute_rrna_score([])
        self.assertEqual(obs, 0.1)

        obs = compute_rrna_score([imd2])
        self.assertEqual(obs, 0.5)

        obs = compute_rrna_score([imd1])
        self.assertEqual(obs, 0.6)

    def test_compute_trna_score(self):
        imd = IntervalMetadata(None)

        obs = compute_trna_score([imd])
        self.assertEqual(obs, 0.1)

        for a in ['Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile',
                  'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']:
            imd.add([(0, 12)], metadata={'type': 'tRNA', 'product': 'tRNA-' + a})

        obs = compute_trna_score([imd])
        self.assertEqual(obs, 0.6)

        for a in ['Ala', 'Arg', 'Asn', 'Asp']:
            imd.add([(0, 12)], metadata={'type': 'tRNA', 'product': 'tRNA-' + a})

        obs = compute_trna_score([imd])
        self.assertEqual(obs, 1)

    def test_compute_gene_score(self):
        seqs = get_data_path('pfam.faa')
        for number, exp in [(1, 0.1), (102, 1)]:
            with NamedTemporaryFile() as faa:
                for i, seq in enumerate(read(seqs, format='fasta')):
                    if i == number:
                        break
                    write(seq, into=faa, format='fasta')
                faa.flush()
                obs = compute_gene_score(faa.name)
                self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
