# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from skbio.util import get_data_path
from skbio.metadata import IntervalMetadata, Feature

from micronota.parsers.cmscan import (_cmscan_to_metadata,
                                      _cmscan_to_generator,
                                      _cmscan_sniffer,
                                      CmscanFormatError)
from wheel.signatures import assertTrue


class CmscanIOTests(TestCase):
    def setUp(self):
        self.file_valid = get_data_path('valid.cmscan')
        self.file_invalidOrientation = get_data_path(
            'invalidOrientation.cmscan')
        self.file_charInPos = get_data_path('charInPos.cmscan')
        self.file_startStopSwop = get_data_path('startStopSwop.cmscan')

        self.subsequences = [{
            'MODEL_NAME': 'LSU_rRNA_bacteria',
            'MODEL_ACCESSION': 'RF02541',
            'SEQUENCE_NAME': 'gi|15829254|ref|NC_002695.1|',
            'SEQUENCE_ACCESSION': '-',
            'TYPE_OF_MODEL': 'cm',
            'MODEL_START_POSITION': '1',
            'MODEL_END_POSITION': '2925',
            'STRAND': '+',
            'TRUNCATED': 'no',
            'PASS': '1',
            'GC_CONTENT': '0.53',
            'BIAS': '45.4',
            'BITSCORE': '2889.8',
            'EVALUE': '0',
            'INC': '!',
            'DESCRIPTION': '-'
        }, {
            'MODEL_NAME': 'SSU_rRNA_microsporidia',
            'MODEL_ACCESSION': 'RF02542',
            'SEQUENCE_NAME': 'gi|15829254|ref|NC_002695.1|',
            'SEQUENCE_ACCESSION': '-',
            'TYPE_OF_MODEL': 'cm',
            'MODEL_START_POSITION': '1',
            'MODEL_END_POSITION': '1312',
            'STRAND': '+',
            'TRUNCATED': 'no',
            'PASS': '1',
            'GC_CONTENT': '0.55',
            'BIAS': '14.5',
            'BITSCORE': '733.3',
            'EVALUE': '5.8e-220',
            'INC': '!',
            'DESCRIPTION': '-'
        }, {
            'MODEL_NAME': 'IS061',
            'MODEL_ACCESSION': 'RF00115',
            'SEQUENCE_NAME': 'gi|secondSeq|',
            'SEQUENCE_ACCESSION': '-',
            'TYPE_OF_MODEL': 'cm',
            'MODEL_START_POSITION': '1',
            'MODEL_END_POSITION': '180',
            'STRAND': '-',
            'TRUNCATED': 'no',
            'PASS': '1',
            'GC_CONTENT': '0.44',
            'BIAS': '0.0',
            'BITSCORE': '232.1',
            'EVALUE': '1.5e-42',
            'INC': '!',
            'DESCRIPTION': '-'
        }, {
            'MODEL_NAME': 'SIB_RNA',
            'MODEL_ACCESSION': 'RF00113',
            'SEQUENCE_NAME': 'gi|secondSeq|',
            'SEQUENCE_ACCESSION': '-',
            'TYPE_OF_MODEL': 'cm',
            'MODEL_START_POSITION': '1',
            'MODEL_END_POSITION': '147',
            'STRAND': '+',
            'TRUNCATED': 'no',
            'PASS': '1',
            'GC_CONTENT': '0.40',
            'BIAS': '0.0',
            'BITSCORE': '136.7',
            'EVALUE': '4e-27',
            'INC': '!',
            'DESCRIPTION': '-'
        }]
        self.intervals = [
            [(4977823, 4980727)],
            [(4831659, 4833190)],
            [(1915102, 1915281)],
            [(3794337, 3794487)]
        ]

        self.seq1 = IntervalMetadata(features={
            Feature(**self.subsequences[0]): self.intervals[0],
            Feature(**self.subsequences[1]): self.intervals[1]
        })
        self.seq2 = IntervalMetadata(features={
            Feature(**self.subsequences[2]): self.intervals[2],
            Feature(**self.subsequences[3]): self.intervals[3]
        })


class ReaderTests(CmscanIOTests):
    def test_cmscan_to_metadata(self):
        t = TestCase()
#       positive controls: will the information from a file parsed to what
#       we expect
        assertTrue(self.seq1.__eq__(_cmscan_to_metadata(
                                                        self.file_valid,
                                                        rec_num=1)))
        assertTrue(self.seq2.__eq__(_cmscan_to_metadata(
                                                        self.file_valid,
                                                        rec_num=2)))

#       negative control: we are parsing the wrong information from the file,
#       check if it is really unequal
        assertTrue(not self.seq1.__eq__(_cmscan_to_metadata(
                                                            self.file_valid,
                                                            rec_num=2)))

#       test if parser raises error about an unknown character as strand
#       identifier
        t.assertRaisesRegex(CmscanFormatError,
                            "Unknown strand character",
                            _cmscan_to_metadata,
                            self.file_invalidOrientation,
                            rec_num=1)

#       test if parser complains about non digit characters in positional
#       arguments
        t.assertRaisesRegex(CmscanFormatError,
                            "must be an integer value for the start position "
                            "of the hit. Here, it is",
                            _cmscan_to_metadata,
                            self.file_charInPos,
                            rec_num=1)

#       test if parser checks for wrong start and stop positions of the hit
#       in the query sequence
        t.assertRaisesRegex(CmscanFormatError,
                            "It might be, that this hit is in fact on the "
                            "reverse strand. Please check strand orientation "
                            "and positions",
                            _cmscan_to_metadata,
                            self.file_startStopSwop,
                            rec_num=1)

    def test_cmscan_to_generator(self):
        assertTrue(list(_cmscan_to_generator(self.file_valid))[0].__eq__(
            self.seq1))
        assertTrue(list(_cmscan_to_generator(self.file_valid))[1].__eq__(
            self.seq2))
        assertTrue(not list(_cmscan_to_generator(self.file_valid))[0].__eq__(
            self.seq2))


class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
                                             'charInPos.cmscan',
                                             'invalidOrientation.cmscan',
                                             'startStopSwop.cmscan',
                                             'valid.cmscan']))
        self.negative_fps = list(map(get_data_path, [
            'blank.sam',
            'uniprot_multi.embl']))

    def test_positive(self):
        for fp in self.positive_fps:
            self.assertEqual(_cmscan_sniffer(fp), (True, {}))
        for fp in self.negative_fps:
            self.assertEqual(_cmscan_sniffer(fp), (False, {}))


if __name__ == '__main__':
    main()
