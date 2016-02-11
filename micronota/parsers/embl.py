r'''

ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt

Classes of entries:

  Class          Definition
  -----------    -----------------------------------------------------------
  CON		 Entry constructed from segment entry sequences; if unannotated,
                 annotation may be drawn from segment entries
  PAT            Patent
  EST            Expressed Sequence Tag
  GSS            Genome Survey Sequence
  HTC            High Thoughput CDNA sequencing
  HTG            High Thoughput Genome sequencing
  MGA            Mass Genome Annotation
  WGS            Whole Genome Shotgun
  TSA            Transcriptome Shotgun Assembly
  STS            Sequence Tagged Site
  STD            Standard (all entries not classified as above)


Taxonomic division
                          Division                 Code
                          -----------------        ----
                          Bacteriophage            PHG
                          Environmental Sample     ENV
                          Fungal                   FUN
                          Human                    HUM
                          Invertebrate             INV
                          Other Mammal             MAM
                          Other Vertebrate         VRT
                          Mus musculus             MUS
                          Plant                    PLN
                          Prokaryote               PRO
                          Other Rodent             ROD
                          Synthetic                SYN
                          Transgenic               TGN
                          Unclassified             UNC
                          Viral                    VRL

structure of an entry

     ID - identification             (begins each entry; 1 per entry)
     AC - accession number           (>=1 per entry)
     PR - project identifier         (0 or 1 per entry)
     DT - date                       (2 per entry)
     DE - description                (>=1 per entry)
     KW - keyword                    (>=1 per entry)
     OS - organism species           (>=1 per entry)
     OC - organism classification    (>=1 per entry)
     OG - organelle                  (0 or 1 per entry)
     RN - reference number           (>=1 per entry)
     RC - reference comment          (>=0 per entry)
     RP - reference positions        (>=1 per entry)
     RX - reference cross-reference  (>=0 per entry)
     RG - reference group            (>=0 per entry)
     RA - reference author(s)        (>=0 per entry)
     RT - reference title            (>=1 per entry)
     RL - reference location         (>=1 per entry)
     DR - database cross-reference   (>=0 per entry)
     CC - comments or notes          (>=0 per entry)
     AH - assembly header            (0 or 1 per entry)
     AS - assembly information       (0 or >=1 per entry)
     FH - feature table header       (2 per entry)
     FT - feature table data         (>=2 per entry)
     XX - spacer line                (many per entry)
     SQ - sequence header            (1 per entry)
     CO - contig/construct line      (0 or >=1 per entry)
     bb - (blanks) sequence data     (>=1 per entry)
     // - termination line           (ends each entry; 1 per entry)

TODO
----
* merge with genbank parsers
  * _parse_genbanks -> _parse_records
* sanity check the length and sequence type with ID line
* parse Reference related lines
* parse CC lines
* parse the DR lines more thoroughly

Examples
--------

Reading UniProt EMBL files
^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> uniprot = [
... r'ID   001R_FRG3G              Reviewed;         256 AA.\n',
... r'AC   Q6GZX4;\N',
... r'DT   28-JUN-2011, integrated into UniProtKB/Swiss-Prot.\n',
... r'DT   19-JUL-2004, sequence version 1.\n',
... r'DT   01-APR-2015, entry version 30.\n',
... r'DE   RecName: Full=Putative transcription factor 001R;\n',
... r'GN   ORFNames=FV3-001R;\n',
... r'OS   Frog virus 3 (isolate Goorha) (FV-3).\n',
... r'OC   Viruses; dsDNA viruses, no RNA stage; Iridoviridae; Ranavirus.\n',
... r'OX   NCBI_TaxID=654924;\n',
... r'OH   NCBI_TaxID=8295; Ambystoma (mole salamanders).\n',
... r'OH   NCBI_TaxID=45438; Rana sylvatica (Wood frog).\n',
... r'RN   [1]\N',
... r'RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].\N',
... r'RX   PubMed=15165820; DOI=10.1016/j.virol.2004.02.019;\n',
... r'RA   Tan W.G., Barkman T.J., Gregory Chinchar V., Essani K.;\n',
... r'RT   "Comparative genomic analyses of frog virus 3, type species of\n',
... r'RT   the genus Ranavirus (family Iridoviridae).";\n',
... r'RL   Virology 323:70-84(2004).\n',
... r'CC   -!- FUNCTION: Transcription activation. {ECO:0000305}.\n',
... r'DR   EMBL; AY548484; AAT09660.1; -; Genomic_DNA.\n',
... r'DR   RefSeq; YP_031579.1; NC_005946.1.\n',
... r'DR   ProteinModelPortal; Q6GZX4; -.\n',
... r'DR   GeneID; 2947773; -.\n',
... r'DR   KEGG; vg:2947773; -.\n',
... r'DR   Proteomes; UP000008770; Genome.\n',
... r'DR   GO; GO:0046782; P:regulation of transcription; IEA:InterPro.\n',
... r'DR   GO; GO:0006351; P:transcription, DNA-templated; IEA:UniProtKB-KW.\n',
... r'DR   InterPro; IPR007031; Poxvirus_VLTF3.\n',
... r'DR   Pfam; PF04947; Pox_VLTF3; 1.\n',
... r'PE   4: Predicted;\n',
... r'KW   Activator; Complete proteome; Reference proteome; Transcription;\n',
... r'KW   Transcription regulation.\n',
... r'FT   CHAIN         1    256       Putative transcription factor 001R.\n',
... r'FT                                /FTId=PRO_0000410512.\n',
... r'FT   COMPBIAS     14     17       Poly-Arg.\n',
... r'SQ   SEQUENCE   256 AA;  29735 MW;  B4840739BF7D4121 CRC64;\N',
... r'     MAFSAEDVLK EYDRRRRMEA LLLSLYYPND RKLLDYKEWS PPRVQVECPK APVEWNNPPS\n',
... r'     EKGLIVGHFS GIKYKGEKAQ ASEVDVNKMC CWVSKFKDAM RRYQGIQTCK IPGKVLSDLD\n',
... r'     AKIKAYNLTV EGVEGFVRYS RVTKQHVAAF LKELRHSKQY ENVNLIHYIL TDKRVDIQHL\n',
... r'     EKDLVKDFKA LVESAHRMRQ GHMINVKYIL YQLLKKHGHG PDGPDILTVK TGSKGVLYDD\n',
... r'     SFRKIYTDLG WKFTPL\N',
... r'//\n']
>>> from skbio import Protein
>>> pro = Protein.read(uniprot)
'''

from skbio.io import create_format
from skbio.sequence import Sequence, DNA, RNA, Protein
from skbio.io.format._base import (
    _line_generator, _get_nth_sequence, _too_many_blanks)


embl = create_format('embl')

# This list is ordered. From EMBL specification
_HEADERS = ['ID',  # identification            (begins each entry; 1 per entry)
            'AC',  # accession number          (>=1 per entry)
            'PR',  # project identifier        (0 or 1 per entry)
            'DT',  # date                      (2 per entry)
            'DE',  # description               (>=1 per entry)
            'KW',  # keyword                   (>=1 per entry)
            'OS',  # organism species          (>=1 per entry)
            'OC',  # organism classification   (>=1 per entry)
            'OG',  # organelle                 (0 or 1 per entry)
            'RN',  # reference number          (>=1 per entry)
            'RC',  # reference comment         (>=0 per entry)
            'RP',  # reference positions       (>=1 per entry)
            'RX',  # reference cross-reference (>=0 per entry)
            'RG',  # reference group           (>=0 per entry)
            'RA',  # reference author(s)       (>=0 per entry)
            'RT',  # reference title           (>=1 per entry)
            'RL',  # reference location        (>=1 per entry)
            'DR',  # database cross-reference  (>=0 per entry)
            'CC',  # comments or notes         (>=0 per entry)
            'AH',  # assembly header           (0 or 1 per entry)
            'AS',  # assembly information      (0 or >=1 per entry)
            'FH',  # feature table header      (2 per entry)
            'FT',  # feature table data        (>=2 per entry)
            'SQ',  # sequence header           (1 per entry)
            'CO',  # contig/construct line     (0 or >=1 per entry)
            'bb']  # (blanks) sequence data    (>=1 per entry)


@embl.sniffer()
def _embl_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}
    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    try:
        _parse_id([line])
    except:
        return False, {}
    return True, {}


@embl.reader(None)
def _embl_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh, _parse_single_embl):
        yield _construct(record, constructor, **kwargs)


@embl.reader(Sequence)
def _embl_to_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    print(record)
    return _construct(record, Protein, **kwargs)


@embl.reader(Protein)
def _embl_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_records(fh, _parse_single_embl), seq_num)
    return _construct(record, Protein, **kwargs)


@embl.reader(DNA)
def _embl_to_DNA(fh, seq_num=1, **kwargs):
    pass


@embl.reader(RNA)
def _embl_to_RNA(fh, seq_num=1, **kwargs):
    pass


def _construct(record, constructor=None, **kwargs):
    seq, md, pmd = record
    if constructor is None:
        unit = md['ID']['unit'].lower()
        if unit == 'bp':
            constructor = DNA
        elif unit == 'aa':
            constructor = Protein
    if constructor == RNA:
        return DNA(
            seq, metadata=md, positional_metadata=pmd, **kwargs).transcribe()
    else:
        return constructor(
            seq, metadata=md, positional_metadata=pmd, **kwargs)


def _parse_records(fh, parser):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith('//'):
            yield parser(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_embl(chunks):
    metadata = {}
    positional_metadata = None
    sequence = ''

    def split(s):
        if s[0].isspace():
            return None, s.strip()
        else:
            return s.split(None, 1)

    for header, data in _yield_section(chunks, split, strip=False):
        if header in ['XX', 'RN', 'RC', 'RP', 'RX', 'RG', 'RA', 'RT', 'RL']:
            continue
        parser = _PARSER_TABLE.get(header)
        if parser:
            data = parser(data)
        if header == 'SQ':
            sequence = data
        else:
            metadata[header] = data

    return sequence, metadata, positional_metadata


def _parse_id(lines):
    '''Parse ID line.'''
    items = lines[0].rstrip('.').split(';')
    res = dict()
    # quality would be "reviewed" or "unreviewed"
    res['id'], res['quality'] = items[0].split()
    res['size'], res['unit'] = items[-1].split()
    return res


def _parse_ac(lines):
    return lines[0].rstrip(';')


def _parse_oc(lines):
    return ' '.join(lines).rstrip('.')


def _parse_ox(lines):
    line = lines[0]
    taxon_id = line.replace('NCBI_TaxID=', '').rstrip(';')
    return taxon_id


def _parse_oh(lines):
    '''Parse the OH (Organism Host) line.

    Optional and appears only in viral entries. It indicates the host
    organism(s) that are susceptible to be infected by a virus.

    Example line:
    'NCBI_TaxID=8295; Ambystoma (mole salamanders).'
    '''
    taxa = dict()
    for line in lines:
        taxon_id, taxon = line.split(';')
        taxon_id = taxon_id.replace('NCBI_TaxID=', '')
        taxon = taxon.strip(' .')
        taxa[taxon_id] = taxon
    return taxa


def _parse_dr(lines):
    '''Parse the DR (Database Cross-reference) line.'''
    res = dict()
    for db_name, db_lines in _yield_section(lines, lambda s: s.split('; ', 1)):
        res[db_name] = db_lines
    return res


def _parse_pe(lines):
    return lines[0].split(':')[0]


def _parse_kw(lines):
    l = ' '.join(lines).rstrip('.')
    return l.split('; ')


def _parse_sq(lines):
    return ''.join(''.join(lines[1:]).split())


def _yield_section(lines, split_header, **kwargs):
    '''Yield the lines with the same header.

    Parameters
    ----------
    split_header : function
        It accepts a string of line and returns the header and
        rest of the needed data as a tuple. If no header exists,
        return ``None`` and the data.

    Notes
    -----
    The following example is a valid section::

    DR   GO GO:000001
    DR   GO GO:000002

    This is also a section::

    SQ   ATGCA ATGCA
         ATGCA ATGCA
    '''
    curr = []
    header, _ = split_header(lines[0])

    for line in _line_generator(lines, **kwargs):
        items = split_header(line)

        # if the header is changed, it is a new section
        if items[0] is not None and items[0] != header:
            if curr:
                yield header, curr
                curr = []
                header = items[0]

        curr.append(items[1].strip())
    # don't forget to return the last section in the file
    if curr:
        yield header, curr


_PARSER_TABLE = {
    'ID': _parse_id,
    'AC': _parse_ac,
    'OC': _parse_oc,
    'OX': _parse_ox,
    'OH': _parse_oh,
    'KW': _parse_kw,
    'DR': _parse_dr,
    'PE': _parse_pe,
    'SQ': _parse_sq}
