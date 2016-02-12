# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util import get_data_path
from skbio import Protein

from micronota.parsers.embl import (
    _embl_sniffer, _embl_to_protein, _embl_to_generator)


class EmblIOTests(TestCase):
    def setUp(self):
        self.multi_fp = get_data_path('uniprot_multi.embl')
        self.single_fp = get_data_path('uniprot_single.embl')
        self.single_exp = (
            'MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGL'
            'IVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNL'
            'TVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALV'
            'ESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL',
            {'AC': 'Q6GZX4',
             'CC': ['-!- FUNCTION: Transcription activation. {ECO:0000305}.',
                    '-----------------------------------------------------------------------',
                    'Copyrighted by the UniProt Consortium, see http://www.uniprot.org/terms',
                    'Distributed under the Creative Commons Attribution-NoDerivs License',
                    '-----------------------------------------------------------------------'],
             'DE': ['RecName: Full=Putative transcription factor 001R;'],
             'DR': {'EMBL': ['AY548484; AAT09660.1; -; Genomic_DNA.'],
                    'GO': ['GO:0046782; P:regulation of viral transcription; IEA:InterPro.',
                           'GO:0006351; P:transcription, DNA-templated; IEA:UniProtKB-KW.'],
                    'GeneID': ['2947773; -.'],
                    'InterPro': ['IPR007031; Poxvirus_VLTF3.'],
                    'KEGG': ['vg:2947773; -.'],
                    'Pfam': ['PF04947; Pox_VLTF3; 1.'],
                    'ProteinModelPortal': ['Q6GZX4; -.'],
                    'Proteomes': ['UP000008770; Genome.'],
                    'RefSeq': ['YP_031579.1; NC_005946.1.']},
             'DT': ['28-JUN-2011, integrated into UniProtKB/Swiss-Prot.',
                    '19-JUL-2004, sequence version 1.',
                    '01-APR-2015, entry version 30.'],
             'FT': ['CHAIN         1    256       Putative transcription factor 001R.',
                    '/FTId=PRO_0000410512.',
                    'COMPBIAS     14     17       Poly-Arg.'],
             'GN': ['ORFNames=FV3-001R;'],
             'ID': {'id': '001R_FRG3G',
                    'quality': 'Reviewed',
                    'size': '256',
                    'unit': 'AA'},
             'KW': ['Activator',
                    'Complete proteome',
                    'Reference proteome',
                    'Transcription',
                    'Transcription regulation'],
             'OC': 'Viruses; dsDNA viruses, no RNA stage; Iridoviridae; Ranavirus',
             'OH': {'45438': 'Rana sylvatica (Wood frog)',
                    '8295': 'Ambystoma (mole salamanders)'},
             'OS': ['Frog virus 3 (isolate Goorha) (FV-3).'],
             'OX': '654924',
             'PE': '4'})


class SnifferTests(TestCase):
    def setUp(self):
        self.fps = list(map(get_data_path, [
            'uniprot_multi.embl',
            'uniprot_single.embl']))

    def test_positive(self):
        for fp in self.fps:
            self.assertEqual(_embl_sniffer(fp), (True, {}))


class ReaderTests(EmblIOTests):
    def test_embl_to_protein(self):
        obs = _embl_to_protein(self.single_fp)
        exp = Protein(self.single_exp[0], self.single_exp[1])
        self.assertEqual(exp, obs)

    def test_embl_to_generator(self):
        exp = Protein(self.single_exp[0], self.single_exp[1])
        for obs in _embl_to_generator(self.multi_fp):
            self.assertEqual(exp, obs)


if __name__ == '__main__':
    main()
