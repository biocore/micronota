from unittest import TestCase, main

from skbio.util import get_data_path
from skbio import Protein, Sequence

from micronota.parsers.sam import (
    _sam_sniffer, _sam_to_protein, _sam_to_generator)


class SamIOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('ecoli_multi.sam')
        self.single_fp = get_data_path('ecoli_single.sam')
        self.single_exp = \
            ('DEQTKRMIRRAILKAVAIPGYQVPFGGREMPMPYGWGTGGIQL'
             'TASVIGESDVLKVIDQGADDTTNAVSIRNFFKRVTGVNTTERT'
             'DDATLIQTRHRIPETPLTEDQIIIFQVPIPEPLRFIEPRETET'
             'RTMHALEEYGVMQVKLYEDIARFGHIATTYAYPVKVNGRYVMD'
             'PSPIPKFDNPKMDMMPALQLFGAGREKRIYAVPPFTHVESLDF'
             'DDHPFTVQQWDEPCAICGSTHSYLDEVVLDDAGNRMFVCSDTD'
             'YCRQQNEAKSQ', {
                 'HD': ['VN:1.5', 'SO:query'],
                 'PG': 'PN:DIAMOND',
                 'mm': 'BlastX',
                 'QNAME': 'WP_000002278.1',
                 'FLAG': 0,
                 'RNAME': 'UniRef100_P16688',
                 'POS': 1,
                 'MAPQ': 255,
                 'CIGAR': '281M',
                 'RNEXT': '*',
                 'PNEXT': 0,
                 'TLEN': 0,
                 'QUAL': '*',
                 'AS': 573,
                 'NM': 3,
                 'ZR': 1477,
                 'ZE': 5.9e-164,
                 'ZI': 98,
                 'ZL': 281,
                 'ZF': 1,
                 'ZS': 1,
                 'MD': '102V117R54S5'
             })


class SnifferTests(SamIOTests):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'ecoli_multi.sam',
            'ecoli_single.sam']))
        self.negative_fps = list(map(get_data_path, [
            'blank.sam']))

    def test_positive(self):
        for fp in self.positive_fps:
            self.assertEqual(_sam_sniffer(fp), (True, {}))

    def test_negative(self):
        for fp in self.negative_fps:
            self.assertEqual(_sam_sniffer(fp), (False, {}))


class ReaderTests(SamIOTests):
    def test_sam_to_protein(self):
        self.maxDiff = None
        obs = _sam_to_protein(self.single_fp)
        exp = Protein(self.single_exp[0],
                      self.single_exp[1])
        self.assertEqual(sorted(obs.metadata.items()),
                         sorted(exp.metadata.items()))

        # FIXME: The equality method in the Sequence object
        # is broken :(
        # self.assertEqual(obs, exp)

    def test_sam_to_generator(self):
        exp = Sequence(self.single_exp[0], self.single_exp[1])
        for obs in _sam_to_generator(self.multi_fp):
            self.assertEqual(sorted(obs.metadata.items()),
                             sorted(exp.metadata.items()))

if __name__ == "__main__":
    main()
