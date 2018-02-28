r'''
GFF3 Parser
=====================

GFF is a standard file format for storing genomic features in a text file.
GFF stands for Generic Feature Format. GFF files are plain text, 9 column,
tab-delimited files [#]_.

The first line of the file is a comment that identifies the format and version.
This is followed by a series of data lines, each one of which corresponds
to an annotation.
The 9 columns of the annotation section are as follows:

    +--------------------------------+----------------------------------+
    | SEQID - ID of the landmark used| (begins each entry; 1 per entry) |
    +--------------------------------+----------------------------------+
    | SOURCE - algorithms used       | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | TYPE - type of the feature     | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | START - start of the feature   | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | END - end of the feature       | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | SCORE - floating point score   | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | STRAND - +/-/./?               | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | PHASE - only for TYPE="CDS"    | 1 per entry                      |
    +--------------------------------+----------------------------------+
    | ATTR - feature attributes (tag)| 1 per entry                      |
    +--------------------------------+----------------------------------+


Column 9 (ATTR) tags have predefined meanings:

ID
Indicates the unique identifier of the feature.
IDs must be unique within the scope of the GFF file.

Name
Display name for the feature. This is the name to be displayed to the user.
Unlike IDs, there is no requirement that the Name be unique within the file.

Alias
A secondary name for the feature.
It is suggested that this tag be used whenever a secondary identifier for
the feature is needed, such as locus names and accession numbers.
Unlike ID, there is no requirement that Alias be unique within the file.

Parent
Indicates the parent of the feature.
A parent ID can be used to group exons into transcripts, transcripts into genes
and so forth. A feature may have multiple parents. Parent can *only* be used
to indicate a partof relationship.

Target
Indicates the target of a nucleotide-to-nucleotide or
protein-to-nucleotide alignment. The format of the value is "target_id start
end [strand]", where strand is optional and may be "+" or "-".
If the target_id contains spaces, they must be escaped as hex escape %20.

Gap
The alignment of the feature to the target if the two are not collinear
(e.g. contain gaps).
The alignment format is taken from the CIGAR format described in the Exonerate
documentation.

Derives_from
Used to disambiguate the relationship between one feature and another when the
relationship is a temporal one rather than a purely structural "part of" one.
This is needed for polycistronic genes.

Note
A free text note.

Dbxref
A database cross reference. See the GFF3 specification for more information.

Ontology_term
A cross reference to an ontology term.

Multiple attributes of the same type are indicated by separating
the values with the comma ","


Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.metadata.IntervalMetadata` objects                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |generator of :mod:`skbio.metadata.IntervalMetadata` objects    |
+------+------+---------------------------------------------------------------+


Reference
---------
.. [#] http://gmod.org/wiki/GFF3
'''

from skbio.io import create_format, FileFormatError
from skbio.metadata import IntervalMetadata, Feature
from skbio.io.format._base import (
    _line_generator, _too_many_blanks)
from skbio.io.format._base import _get_nth_sequence as _get_nth_record


class GFFFormatError(FileFormatError):
    pass

gff = create_format('gff')

# Annotation headers
_ANNOTATION_HEADERS = [
          'SEQID',
          'SOURCE',
          'TYPE',
          'SCORE',
          'STRAND',
          'PHASE',
          'ATTR'
          ]

_GFF_HEADERS = [
          'SEQID',
          'SOURCE',
          'TYPE',
          'START',
          'END',
          'SCORE',
          'STRAND',
          'PHASE',
          'ATTR'
          ]

# Attribute headers
_ATTR_HEADERS = [
        'ID',
        'Name',
        'Alias',
        'Parent',
        'Target',
        'Gap',
        'Derives_from',
        'Note',
        'Dbxref',
        'Ontology_term'
        ]


@gff.sniffer()
def _gff_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    if line.startswith('##gff-version 3'):
        return True, {}
    else:
        return False, {}


def _is_float(input):
    try:
        float(input)
    except ValueError:
        return False
    return True


def _parse_attr(s):
    _field = s.split(';')
    _attr = ()
    for f in _field:
        try:
            _type, _val = f.split('=')
        except ValueError:
            _type, _val = 0, 0

        _attr = _attr + (_type, _val)
    return _attr


def _parse_required(s):
    if s.isdigit():
        return int(s)
    elif _is_float(s):
        return float(s)
    else:
        return s


def _construct(record, constructor=None, **kwargs):
    if constructor is None:
        constructor = IntervalMetadata
    if constructor is IntervalMetadata:
        return IntervalMetadata(features=record)


@gff.reader(None)
def _gff_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield _construct(record, constructor, **kwargs)


@gff.reader(IntervalMetadata)
def _gff_to_metadata(fh, rec_num=1, **kwargs):
    record = _get_nth_record(_parse_records(fh), rec_num)
    return _construct(record, IntervalMetadata, **kwargs)


def _parse_records(fh):
    seqID = False
    features = {}
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if not line.startswith('#'):
            tabs = line.split('\t')

            # extract annotations
            _headers = tabs[:3] + tabs[5:8]

            # extract attributes
            _attributes = tabs[8]

            # extract intervals
            _intervals = tabs[3:5]

            intervals = tuple(map(_parse_required, _intervals))

            attr = _parse_attr(_attributes)
            _annot = list(map(_parse_required, _headers))
            _annot.append(attr)
            annotation = Feature(dict(zip(_ANNOTATION_HEADERS, _annot)))

            if _annot[0] == seqID or seqID is False:
                seqID = _annot[0]
                features[annotation] = intervals
            else:
                seqID = False
                yield features
                features = {}
                features[annotation] = intervals

    yield features


def _attr_to_list(attr_list):
    _tags = []
    if any(isinstance(el, tuple) for el in attr_list):
        for _attr in attr_list:
            _tags.append('='.join([_attr[0], _attr[1]]))
        return ';'.join(_tags)
    else:
        return '='.join([attr_list[0], attr_list[1]])


@gff.writer(IntervalMetadata)
def _IntervalMetadata_to_gff(obj, fh):
    # write file header
    fh.write('%s\n' % '##gff-version 3')

    for features, interval in obj.features.items():
        tab = [None] * 9

        if len(features) is not 7:
            raise GFFFormatError(
                "``IntervalMetadata`` can only be written in GFF format if all"
                " annotation columns are found.")
        if len(interval) is not 2:
            raise GFFFormatError(
                "``IntervalMetadata`` can only be written in GFF format if "
                " `START` and `END` fields are provided.")
        if not all(_annot in set(features) for _annot in _ANNOTATION_HEADERS):
            raise GFFFormatError(
                "GFF format requires header names to match pre-defined set: %s"
                % ', '.join(_ANNOTATION_HEADERS))

        for i, annot in enumerate(_GFF_HEADERS):
            if annot is 'START':
                tab[i] = interval[0]
            elif annot is 'END':
                tab[i] = interval[1]
            elif annot is 'ATTR':
                tab[i] = _attr_to_list(features[annot])
            else:
                tab[i] = features[annot]

        fh.write('\t'.join(map(str, tab)) + '\n')
