r'''
TIGRFAM
=======

.. currentmodule:: micronota.db.tigrfam

This module create reference database from TIGRFAM [#]_ for micronota usage.
TIGRFAM may be used with HMMER3. The HMM libraries should be usable in
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
       ``subfamily`` model, which in turn outranks a ``domain`` model [#]_.
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


* TIGRFAMS_GO_LINK

  Gene Ontology (GO) term assignments.
  Models may have  0, 1, or several GO assignments.


Functions
---------

.. autosummary::
   :toctree: generated/

   prepare_db
   prepare_metadata

Reference
---------
.. [#] http://www.ncbi.nlm.nih.gov/pubmed/12520025
.. [#] http://www.ncbi.nlm.nih.gov/pubmed/23197656

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
from os.path import join, basename
from tempfile import mkdtemp
from sqlite3 import connect
from logging import getLogger

from ..bfillings.hmmer import run_hmmpress

from ..util import _overwrite, _download


def prepare_db(out_d, downloaded, prefix='tigrfam_v15.0', force=False,
               hmm='ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz',
               metadata='ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_INFO.tar.gz'):
    '''Download and prepare TIGRFAM database.

    Parameters
    ----------
    out_d : str
        The directory of output files.
    downloaded : str
        The directory of downloaded files.
    prefix : str
        The file name (without extensions) of the output files.
    force : boolean
        Whether to overwrite existing files
    hmm : str
        The file name of hmm models
    metadata : str
        The file name of the metadata for the hmm models
    '''
    logger = getLogger(__name__)
    logger.info('Preparing %s database' % prefix)

    hmm_fp = join(out_d, '%s.hmm' % prefix)
    hmm_raw = join(downloaded, basename(hmm))
    metadata_fp = join(out_d, '%s.db' % prefix)
    metadata_raw = join(downloaded, basename(metadata))
    metadata_dir = mkdtemp()
    try:
        # fetch metadata file
        _download(metadata, metadata_raw, overwrite=force)
        # fetch HMM model file
        _download(hmm, hmm_raw, overwrite=force)
    except FileExistsError:
        pass

    with tarfile.open(metadata_raw) as tar:
        tar.extractall(metadata_dir)
    prepare_metadata(metadata_dir, metadata_fp, force)
    shutil.rmtree(metadata_dir)

    # gunzip and move the file
    with gzip.open(hmm_raw, 'rb') as i_f, open(hmm_fp, 'wb') as o_f:
        shutil.copyfileobj(i_f, o_f)

    # don't forget to compress the hmm file
    run_hmmpress(hmm_fp)


def prepare_metadata(in_d, out_fp, overwrite=True):
    '''Compile the metadata into sqlite3 database.

    Parameters
    ----------
    in_d : str
        The input directory containing XXX.INFO files.
    out_fp : str
        The output file path of sqlite3 database.
    overwrite : boolean
        Whether to overwrite if the ``out_fp` exists.

    Returns
    -------
    int
        The number of records processed.

    Notes
    -----
    The schema of the database file contains one table named `tigrfam` that
    has following columns:

    1. ``ac``. TEXT. TIGRFAM accession.

    2. ``key``. TEXT. 2-letter key described in the table above in this module.

    3. ``val``. BLOB. The value of the key.

    4. ``transfer``. INTEGER. Used as the boolean. ``1`` means the ``val``
       should be transferred to the query sequences as its annotation;
       ``0`` means not.

    The table in the database file will be dropped and re-created if
    the function is re-run.
    '''
    logger = getLogger(__name__)
    logger.info('Preparing metadata db for TIGRFAM')

    n = 0
    _overwrite(out_fp, overwrite)
    with connect(out_fp) as conn:
        table_name = 'metadata'
        conn.execute('''CREATE TABLE IF NOT EXISTS {t} (
                            ac       TEXT    NOT NULL,
                            key      TEXT    NOT NULL,
                            val      BLOB    NOT NULL,
                            transfer BOOLEAN NOT NULL,
                        CHECK (transfer IN (0, 1)))'''.format(t=table_name))

        for f in os.listdir(in_d):
            if f.startswith('.') or not f.endswith('.INFO'):
                continue
            n += 1
            tigrfam_id = f.split('.', 1)[0]
            insert = '''INSERT INTO {t} (ac, key, val, transfer)
                        VALUES (?,?,?,?)'''.format(t=table_name)
            for i, j, k in _read_info(join(in_d, f)):
                conn.execute(insert, (tigrfam_id, i, j, k))
        # don't forget to index the column to speed up query
        conn.execute('CREATE INDEX ac ON {t} (ac);'.format(t=table_name))
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
        key, val, int of 0 or 1.
    '''
    # the error param is for the non-utf8 symbols
    with open(fn, errors='backslashreplace') as f:
        for line in f:
            line = line.strip()
            key = line[:2]
            val = line[2:].strip()
            if key in ['TC', 'NC']:
                g, d = [float(i) for i in val.split()]
                yield '%s_global' % key, g, 0
                yield '%s_domain' % key, d, 0
            elif key == 'EC':
                for n in val.split():
                    yield key, n, 1
            elif key == 'RM':
                s = 'PMID:'
                if val.startswith(s):
                    val = val.replace(s, '').strip()
            elif key in ['TP', 'AC', 'RN', 'RT', 'RA', 'RL']:
                continue
            else:
                yield key, val, 1
