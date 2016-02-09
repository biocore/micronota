r'''
TIGRFAM
=======

.. currentmodule:: micronota.db.tigrfam

TIGRFAMs [1]_ may be used with HMMER3. The HMM libraries should be usable in
the same way as Pfam libraries starting from Pfam release 24.0.


Release
-------
:version: 15.0
:date:    9/17/2014
:link:    ftp://ftp.tigr.org/pub/data/TIGRFAMs/RELEASE_NOTE_15.0


Files
-----
* TIGRFAMs_15.0_HMM.LIB.gz

standard TIGRFAMs HMMs, ASCII. The highest accession number of the current
release, 15.0, is TIGR04571.

* TIGRFAMs_15.0_INFO.tar.gz

metadata for each TIGRFAMs. Each family file contains:

   === ========================================================================
   Tag Description
   === ========================================================================
   ID  Identification: One word less than 16 characters
   AC  Accession number. TIGR accession numbers take the form TIGRxxxxx where
       x is a digit
   DE  description of the HMM
   AU  Author. Person(s) responsible for alignment in the format, eg. Fish N
   TC  Trusted cutoffs: global value, then domain value
   NC  Noise cutoffs: global value, then domain value
   CC  Comment. Comment lines may repeat. This is free text
   RN  Reference number
   RM  Reference Medline now holds mostly (but not quite always) PUBMED ids
       (PMID)
   RT  Reference title
   RA  Reference author line
   RL  Reference location (journal name, volume, pages, year).
   DR  Database references. Some differences in databases cited and in format
       vs. Pfam may occur.
   === ========================================================================

Some tags were introduced to support definition of TIGRFAMs or
their use in annotation:

   === ========================================================================
   Tag Description
   === ========================================================================
   AL  Alignment method of seed
   IT  The "isology type", or homology type. "equivalog" models can be used
       for automatic annotation of protein name, prokaryotic gene symbol, EC
       (enzyme commission) number, and role category.
       An ``equivalog`` model assigns more specific annotations than a
       ``subfamily`` model, which in turn outranks a ``domain`` model [2]_.
   GS  Gene symbol - can be applied automatically for prokaryotic
       sequences.
   EC  Enzyme Commission number. In the format  6.1.1.7  without the EC.
       This field may contain more than one EC number, with a single space as
       the separator.
   EN  Expanded name - A fuller or informative alternate version of the
       definition line. This compares (or is identical) to the more terse DE
       definition that is used for automated annotation.
   TP  Always "TIGRFAMs", the identifier of this database
   === ========================================================================

Proteins that score above the trusted cutoffs are believed to reside within
the family and those falling below the noise cutoffs are believed to reside
outside the family.  The margin of error with respect to presence or absence
of a protein within a TIGRFAMs family is represented by the score range
between noise and trusted cutoffs.

Because the number of completed and nearly completed genomes
has now entered the thousands, protein families are becoming
very large, and exceptions may be found in certain equivalog
families.  Exceptions usually represent neofunctionalizations
that arise within an equivalog family, although some may
represent paralogs that are minimally derived since their
branching from the equivalog familiy (i.e. short branch length).

TIGRFAMs now includes 57 models of type "exception", a type
of model that overrules annotation from an equivalog model,
either to give more specific information or to correct annotation
for a neofunctionalized subgroup.

TIGRFAMS_GO_LINK
^^^^^^^^^^^^^^^^
Gene Ontology (GO) term assignments.
Models may have  0, 1, or several GO assignments.


Reference
---------
.. [1] http://www.ncbi.nlm.nih.gov/pubmed/12520025
.. [2] http://www.ncbi.nlm.nih.gov/pubmed/23197656

'''

# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import gzip
import tarfile
from tempfile import mkdtemp
from os.path import join, exists
from ftplib import FTP
from sqlite3 import connect

from ..bfillings.hmmer import hmmpress_hmm


def prepare_db(out_d, prefix='tigrfam_v15.0', force=False,
               server='ftp.tigr.org',
               path='pub/data/TIGRFAMs',
               hmm='TIGRFAMs_15.0_HMM.LIB.gz',
               metadata='TIGRFAMs_15.0_INFO.tar.gz'):
    '''Download and prepare TIGRFAM database.

    Parameters
    ----------
    out_d : str
        The directory of output files.
    prefix : str
        The file name (without extensions) of the output files.
    force : boolean
        Whether to overwrite existing files
    server : str
        FTP server to download TIGRFAM files
    path : str
        The path on FTP to find files to download
    hmm : str
        The file name of hmm models
    metadata : str
        The file name of the metadata for the hmm models
    '''

    hmm_out = join(out_d, '%s.hmm' % prefix)
    if exists(hmm_out) and not force:
        print('Database %s already exists. Skip it.' % prefix)
        return

    metadata_out = join(out_d, '%s.db' % prefix)

    try:
        temp_dir = mkdtemp()
        hmm_tmp = join(temp_dir, hmm)
        metadata_tmp = join(temp_dir, metadata)
        with FTP(server) as ftp:
            ftp.login()
            ftp.cwd(path)

            # fetch metadata file
            ftp.retrbinary('RETR %s' % metadata,
                           open(metadata_tmp, 'wb').write)
            with tarfile.open(metadata_tmp) as tar:
                tar.extractall(temp_dir)
                prepare_metadata(temp_dir, metadata_out)

            # fetch HMM model file
            ftp.retrbinary('RETR %s' % hmm,
                           open(hmm_tmp, 'wb').write)
            # gunzip and move the file
            with gzip.open(hmm_tmp, 'rb') as i_f, open(hmm_out, 'wb') as o_f:
                shutil.copyfileobj(i_f, o_f)
            # don't forget to compress the hmm model file
            hmmpress_hmm(hmm_out)
    finally:
        shutil.rmtree(temp_dir)


def prepare_metadata(in_d, fp):
    '''Compile the metadata into sqlite3 database.

    Parameters
    ----------
    in_d : str
        The input directory containing XXX.INFO files.
    fp : str
        The output file path of sqlite3 database.

    Returns
    -------
    int
        The number of records processed.

    Notes
    -----
    The schema of the database file contains one table named `tigrfam` that
    has following columns:

    1. ``id``. TEXT. TIGRFAM ID.

    2. ``tag``. TEXT. 2-letter tag described in the table above in this module.

    3. ``value``. BLOB. The value of the tag.

    4. ``transfer``. INTEGER. Used as the boolean. ``1`` means the ``value``
       should be transferred to the query sequences as its annotation;
       ``0`` means not.

    The table in the database file will be dropped and re-created if
    the function is re-run.
    '''
    n = 0
    with connect(fp) as conn:
        conn.execute("DROP TABLE IF EXISTS tigrfam")
        conn.execute('''CREATE TABLE tigrfam (
                            id       TEXT    NOT NULL,
                            tag      TEXT    NOT NULL,
                            value    BLOB    NOT NULL,
                            transfer BOOLEAN NOT NULL,
                        CHECK (transfer IN (0,1)));''')

        for f in os.listdir(in_d):
            if not f.endswith('.INFO'):
                continue
            n += 1
            tigrfam_id = f.split('.', 1)[0]
            for x in _read_info(join(in_d, f)):
                conn.execute('''INSERT INTO tigrfam (id, tag, value, transfer)
                                    VALUES ("%s",?,?,?);''' % tigrfam_id,
                             x)
        # don't forget to index the column to speed up query
        conn.execute('CREATE INDEX id ON tigrfam (id);')
        conn.commit()
    return n


def _read_info(fn):
    '''Parse the .INFO file.

    Parameters
    ----------
    fn : str
        file path

    Yields
    ------
    tuple
        tag, value, int of 0 or 1.
    '''
    # the error param is for the non-utf8 symbols
    with open(fn, errors='backslashreplace') as f:
        for line in f:
            line = line.strip()
            tag = line[:2]
            value = line[2:].strip()
            if tag in ['TC', 'NC']:
                g, d = [float(i) for i in value.split()]
                yield '%s_global' % tag, g, 0
                yield '%s_domain' % tag, d, 0
            elif tag == 'EC':
                for n in value.split():
                    yield tag, n, 1
            elif tag == 'RM':
                s = 'PMID:'
                if value.startswith(s):
                    value = value.replace(s, '').strip()
            elif tag in ['TP', 'AC', 'RN', 'RT', 'RA', 'RL']:
                continue
            else:
                yield tag, value, 1
