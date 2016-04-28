r'''
CMSCAN Parser
=============

The program cmscan of the Infernal suite [1]_ is a tools to screen a nucleotide
sequence for homology to known functional RNA families. Input is a set of fasta
sequences and a set of family models. Output is a set of potential hits of sub-
sequences to one of the family models.

Output of cmscan comes in two flavors. STDOUT provides complete hit
reports, but is hard to be parsed. If the flag --tbloutput FILENAME is used,
FILENAME "save parseable table of hits to file". (Alignment between query
sequence and family model are **not** included.) This python module is for
parsing the tbloutput format.

This example lists the first hits between the chromosome of E.Coli and Rfam 12
using CMscan:

.. code-block:: none
#target name           accession query name                   accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
#--------------------- --------- ---------------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
LSU_rRNA_bacteria      RF02541   gi|15829254|ref|NC_002695.1| -          cm        1     2925  4977823  4980727      +    no    1 0.53  45.4 2889.8         0 !   -
LSU_rRNA_bacteria      RF02541   gi|15829254|ref|NC_002695.1| -          cm        1     2925  5018849  5021753      +    no    1 0.53  45.4 2889.8         0 !   -
LSU_rRNA_bacteria      RF02541   gi|15829254|ref|NC_002695.1| -          cm        1     2925  4737148  4740052      +    no    1 0.53  45.4 2889.5         0 !   -
LSU_rRNA_bacteria      RF02541   gi|15829254|ref|NC_002695.1| -          cm        1     2925  4833642  4836546      +    no    1 0.53  45.4 2889.5         0 !   -
LSU_rRNA_bacteria      RF02541   gi|15829254|ref|NC_002695.1| -          cm        1     2925   229090   231992      +    no    1 0.53  43.8 2888.4         0 !   -
LSU_rRNA_archaea       RF02540   gi|15829254|ref|NC_002695.1| -          cm        1     2990  4977822  4980726      +    no    1 0.53  45.8 1848.7         0 !   -


Lines starting with \# are considered as comments and should not be trusted at
all, i.e. the column names in the first line are just for human readability and
the dashes in the second line do not specify column width, as one could wrongly
assume. For all other lines, columns are separated by one (or several) spaces,
since spaces are not allowed to occure within a column, except the last one
(personal communication with Sean R Eddy). Thus, we need to know that there are
18 columns (from the Infernal Userguide, section 6):

(Note that *target* and *query* is used differently between CMsearch and
CMscan. For CMscan, the target is the family and the query is the sequence. It
is the exact opposite for CMsearch.)

1. target name
   The name of the target sequence or profile.
2. accession
   The accession of the target sequence or profile, or '-' if none.
3. query name
   The name of the query sequence or profile.
4. accession
   The accession of the query sequence or profile, or '-' if none.
5. mdl (model)
   Which type of model was used to compute the final score. Either 'cm' or
   'hmm'. A CM is used to compute the final hit scores unless the model has
   zero basepairs or the --hmmonly option is used, in which case a HMM will be
   used.
6. mdl from (model coord)
   The start of the alignment of this hit with respect to the profile (CM or
   HMM), numbered 1..N for a profile of N consensus positions.
7. mdl to (model coord)
   The end of the alignment of this hit with respect to the profile (CM or
   HMM), numbered 1..N for a profile of N consensus positions.
8. seq from (ali coord)
   The start of the alignment of this hit with respect to the sequence,
   numbered 1..L for a sequence of L residues.
9. seq to (ali coord)
   The end of the alignment of this hit with respect to the sequence, numbered
   1..L for a sequence of L residues.
10. strand
   The strand on which the hit occurs on the sequence. '+' if the hit is on the
   top (Watson) strand, '-' if the hit is on the bottom (Crick) strand. If on
   the top strand, the "seq from" value will be less than or equal to the
   "seq to" value, else it will be greater than or equal to it.
11. trunc
   Indicates if this is predicted to be a truncated CM hit or not. This will be
   "no" if it is a CM hit that is not predicted to be truncated by the end of
   the sequence, "5' " or "3' " if the hit is predicted to have one or more 5'
   or 3' residues missing due to a artificial truncation of the sequence, or
   "5'&3'" if the hit is predicted to have one or more 5' residues missing and
   one or more 3' residues missing. If the hit is an HMM hit, this will always
   be '-'.
12. pass
   Indicates what "pass" of the pipeline the hit was detected on. This is
   probably only useful for testing and debugging. Non-truncated hits are found
   on the first pass, truncated hits are found on successive passes.
13. gc
   Fraction of G and C nucleotides in the hit.
14. bias
   The biased-composition correction: the bit score difference contributed by
   the null3 model for CM hits, or the null2 model for HMM hits. High bias
   scores may be a red flag for a false positive. It is difficult to correct
   for all possible ways in which a nonrandom but nonhomologous biological
   sequences can appear to be similar, such as short-period tandem repeats, so
   there are cases where the bias correction is not strong enough (creating
   false positives).
15. score
   The score (in bits) for this target/query comparison. It includes the
   biased-composition correction (the "null3" model for CM hits, or the "null2"
   model for HMM hits).
16. E-value
   The expectation value (statistical significance) of the target. This is a
   per query E-value; i.e. calculated as the expected number of false positives
   achieving this comparison's score for a single query against the search
   space Z. For cmsearch Z is defined as the total number of nucleotides in the
   target dataset multiplied by 2 because both strands are searched. For cmscan
   Z is the total number of nucleotides in the query sequence multiplied by 2
   because both strands are searched and multiplied by the number of models in
   the target database. If you search with multiple queries and if you want to
   control the overall false positive rate of that search rather than the false
   positive rate per query, you will want to multiply this per-query E-value by
   how many queries you're doing.
17. inc
   Indicates whether or not this hit achieves the inclusion threshold: '!' if
   it does, '?' if it does not (and rather only achieves the reporting
   threshold). By default, the inclusion threshold is an E-value of 0.01 and
   the reporting threshold is an E-value of 10.0, but these can be changed with
   command line options as described in the manual pages.
18. description of target
   The remainder of the line is the target's description line, as free text.


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
.. [1] Eric P. Nawrocki and Sean R. Eddy, "Infernal 1.1: 100-fold faster RNA
       homology searches",  Bioinformatics 2013,
       doi: 10.1093/bioinformatics/btt509

'''

from skbio.io import create_format, FileFormatError
from skbio.metadata import IntervalMetadata, Feature
from skbio.io.format._base import (_line_generator, _too_many_blanks)
from skbio.io.format._base import _get_nth_sequence as _get_nth_record

cmscan = create_format('cmscan')

# column headers
_COLUMNS = {
   'MODEL_NAME':              {'position':  0, 'width': 22, 'textalignment': 'l', 'caption': 'target name'},
   'MODEL_ACCESSION':         {'position':  1, 'width':  9, 'textalignment': 'l', 'caption': 'accession'},
   'SEQUENCE_NAME':           {'position':  2, 'width': 28, 'textalignment': 'l', 'caption': 'query name'},
   'SEQUENCE_ACCESSION':      {'position':  3, 'width':  9, 'textalignment': 'l', 'caption': 'accession'},
   'TYPE_OF_MODEL':           {'position':  4, 'width':  3, 'textalignment': 'r', 'caption': 'mdl'},
   'MODEL_START_POSITION':    {'position':  5, 'width':  8, 'textalignment': 'r', 'caption': 'mdl from'},
   'MODEL_END_POSITION':      {'position':  6, 'width':  8, 'textalignment': 'r', 'caption': 'mdl to'},
   'SEQUENCE_START_POSITION': {'position':  7, 'width':  8, 'textalignment': 'r', 'caption': 'seq from'},
   'SEQUENCE_END_POSITION':   {'position':  8, 'width':  8, 'textalignment': 'r', 'caption': 'seq to'},
   'STRAND':                  {'position':  9, 'width':  6, 'textalignment': 'r', 'caption': 'strand'},
   'TRUNCATED':               {'position': 10, 'width':  5, 'textalignment': 'r', 'caption': 'trunc'},
   'PASS':                    {'position': 11, 'width':  4, 'textalignment': 'r', 'caption': 'pass'},
   'GC_CONTENT':              {'position': 12, 'width':  4, 'textalignment': 'r', 'caption': 'gc'},
   'BIAS':                    {'position': 13, 'width':  5, 'textalignment': 'r', 'caption': 'bias'},
   'BITSCORE':                {'position': 14, 'width':  6, 'textalignment': 'r', 'caption': 'score'},
   'EVALUE':                  {'position': 15, 'width':  9, 'textalignment': 'r', 'caption': 'E-value'},
   'INC':                     {'position': 16, 'width':  3, 'textalignment': 'l', 'caption': 'inc'},
   'DESCRIPTION':             {'position': 17, 'width': 21, 'textalignment': 'l', 'caption': 'description of target'}
}
_orderedKeys = sorted(_COLUMNS, key=lambda x: _COLUMNS[x]['position'])


class CmscanFormatError(FileFormatError):
    pass


def _construct(record, constructor=None, **kwargs):
    if constructor is None:
        constructor = IntervalMetadata
    if constructor is IntervalMetadata:
        return IntervalMetadata(features=record)


@cmscan.sniffer()
def _cmscan_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    if line.startswith('#target name'):
        return True, {}
    else:
        return False, {}


@cmscan.reader(None)
def _cmscan_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_records(fh):
        yield _construct(record, constructor, **kwargs)


@cmscan.reader(IntervalMetadata)
def _cmscan_to_metadata(fh, rec_num=1, **kwargs):
    record = _get_nth_record(_parse_records(fh), rec_num)
    return _construct(record, IntervalMetadata, **kwargs)


def _parse_records(fh):
    """ parses an tabular output file, generated by either cmsearch or cmscan
    of the Infernal package, to a collection of IntervalMetadata objects

    Parameters
    ----------
    fh : file handle 
         A file handle to the file that should be parsed
    
    Returns
    -------
    A collection of IntervalMetadata objects
    """
    program = 'CMSCAN'
    currentsequence = False
    annotations = {}
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if not line.startswith('#'):
            attributes = {}
            fields = line.split()  # split at multiple occurrences of spaces

            # checking data for start end end position of the hit in the query
            # sequence.
            hitStart = -1
            if fields[_COLUMNS['SEQUENCE_START_POSITION']['position']].isdigit():
                hitStart = int(
                   fields[_COLUMNS['SEQUENCE_START_POSITION']['position']]
                )
            else:
                raise CmscanFormatError("%s %i %s '%s'." % (
                       "Column",
                       _COLUMNS['SEQUENCE_START_POSITION']['position'],
                       "must be an integer value for the start position of the"
                       " hit. Here, it is",
                       fields[_COLUMNS['SEQUENCE_START_POSITION']['position']])
                      )

            hitEnd = -1
            if fields[_COLUMNS['SEQUENCE_END_POSITION']['position']].isdigit():
                hitEnd = int(
                   fields[_COLUMNS['SEQUENCE_END_POSITION']['position']]
                )
            else:
                raise CmscanFormatError("%s %i %s '%s'." % (
                        "Column",
                        _COLUMNS['SEQUENCE_END_POSITION']['position'],
                        "must be an integer value for the end position of the "
                        "hit. Here, it is",
                        fields[_COLUMNS['SEQUENCE_END_POSITION']['position']],
                        ".")
                      )

            hitOrientation = fields[_COLUMNS['STRAND']['position']]
            if hitOrientation == "+":
                if hitStart > hitEnd:
                    raise CmscanFormatError('%s %i %s %i %s' % (
                          "On the forward strand (+), start position of a hit "
                          "must always be smaller than its end position."
                          " This is not true for the hit between",
                          hitStart,
                          "and",
                          hitEnd,
                          ". It might be, that this hit is in fact on the "
                          "reverse strand. Please check strand orientation and"
                          " positions."))
            elif hitOrientation == "-":
                if hitStart < hitEnd:
                    raise CmscanFormatError('%s %i %s %i %s' % (
                          "On the reverse strand (-), start position of a hit "
                          "must always be larger than its end position."
                          " This is not true for the hit between",
                          hitStart,
                          "and",
                          hitEnd,
                          ". It might be, that this hit is in fact on the "
                          "forward strand. Please check strand orientation and"
                          " positions."))
                else:
                    hitStart, hitEnd = hitEnd, hitStart  # swap orientation
            else:
                raise CmscanFormatError(
                    "%s '%s' %s %i. %s" %
                    ("Unknown strand character",
                     hitOrientation, "in column",
                     _COLUMNS['STRAND']['position'],
                     "Valid characters are '+' for the forward strand and "
                     "'-' for the reverse strand."))

            # since Infernal want the user to be aware of differences between
            # CMsearch and CMscan the information about model and query are at
            # different columns. Here, we contradict this design and store the
            # information always at the same key
            if program is 'CMSEARCH':
                attributes['SEQUENCE_NAME'] = fields[0]
                attributes['SEQUENCE_ACCESSION'] = fields[1]
                attributes['MODEL_NAME'] = fields[2]
                attributes['MODEL_ACCESSION'] = fields[3]
            elif program is 'CMSCAN':
                attributes['SEQUENCE_NAME'] = fields[2]
                attributes['SEQUENCE_ACCESSION'] = fields[3]
                attributes['MODEL_NAME'] = fields[0]
                attributes['MODEL_ACCESSION'] = fields[1]
            else:
                raise CmscanFormatError("Argument 'program' must be either "
                                        "'CMSEARCH' or 'CMSCAN'!")

            # iterate through all keys that have not already be handled above
            for key in _orderedKeys:
                if key in ['SEQUENCE_START_POSITION',
                           'SEQUENCE_END_POSITION',
                           'SEQUENCE_NAME',
                           'SEQUENCE_ACCESSION',
                           'MODEL_NAME',
                           'MODEL_ACCESSION']:
                    continue
                attributes[key] = fields[_COLUMNS[key]['position']]

            # cmscan works on multiple sequence in one FASTA file. We want to
            # yield a separate object for each sequence, thus we create a new
            # one whenever the ID changes
            if (currentsequence != attributes['SEQUENCE_NAME'] and
                    currentsequence is not False):
                yield annotations
                annotations = {}

            # a metadata interval is made out of 'frozen dictionary' aka
            # Features object, where the its hashable value constitutes the
            # key and the value is a list of tuples, aka intervals
            annotations[Feature(attributes)] = [(hitStart, hitEnd)]
            # store current sequence id for the next iteration
            currentsequence = attributes['SEQUENCE_NAME']

    yield annotations


def _writeheader(fh):
    """ Print a header for cmscan hits. It is similar to 'CMscan --tblout' but
        not identical, since we do not check for longest occurring names, but
        use fixed field width. """

    # write file header
    lineCaption = '#'
    lineDelimeter = '#'
    for key in _orderedKeys:
        length = _COLUMNS[key]['width']
        if _COLUMNS[key]['position'] == 0:
            length = length-1
        lineCaption = lineCaption + _printField(
               _COLUMNS[key], _COLUMNS[key]['caption'], length)
        lineDelimeter = lineDelimeter + ('-' * length)
        if _COLUMNS[key]['position']+1 != len(_COLUMNS):
            lineCaption = lineCaption + " "
            lineDelimeter = lineDelimeter + " "
    fh.write(lineCaption + "\n")
    fh.write(lineDelimeter + "\n")


@cmscan.writer(IntervalMetadata)
def _IntervalMetadata_to_cmscan(obj, fh):
    """ Prints IntervalMetadata object in a format similar to CMscan --tblout.
        It is not identical, since ordering is arbitrary here. Furthermore,
        field width might vary. """
    # write data
    for hit, interval in obj.features.items():
        if len(interval) != 1:
            raise CmscanFormatError(
             "The cmscan format allows only one interval per IntervalMetadata "
             "object. Cannot continue writing output.")
        lineData = ''
        for key in _orderedKeys:
            if key == 'SEQUENCE_START_POSITION':
                if hit['STRAND'] == '+':
                    value = str(interval[0][0])
                else:
                    value = str(interval[0][1])
            elif key == 'SEQUENCE_END_POSITION':
                if hit['STRAND'] == '+':
                    value = str(interval[0][1])
                else:
                    value = str(interval[0][0])
            else:
                value = hit[key]
            lineData = lineData + _print_field(_COLUMNS[key], value)
            if _COLUMNS[key]['position']+1 != len(_COLUMNS):
                lineData = lineData + " "
        fh.write(lineData + "\n")


def _print_field(column, text, width=-1):
    """ Helper function to print a string in a column of specific width either
    aligned to the left or the right.

    Parameters
    ----------
    column : dict
            A dictionary holding information about the column to be printed,
            c.f. the _COLUMN attribute. The dict must hold at least the keys 
            'width' and 'textalignment', where width is the column
            width in terms of number of characters, and alignment must be either
            the character "l" or "r" to indicate left or right text alignment 
            within the column 
    text :  string
            the text to be printed
    width : int
            manually override the width value of the 'column' dictionary

    Return
    ------
    A string with leading (r) or trailing (l) spaces and the text if not longer
    than specified width. Otherwise text will be truncated to width."""
    if width < 0:
        width = column['width']
    space = ' ' * (width-len(text))
    if field['textalignment'] == 'l':
        return "%s%s" % (text, space)
    else:
        return "%s%s" % (space, text)
