# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger
from os.path import exists, join, basename, splitext
from subprocess import CalledProcessError

from dumpling import check_range, Dumpling, ArgmntParam, OptionParam, Parameters


params = [
    OptionParam('-searchWL',
                action=check_range(6, 9),
                help='Length of search window used to discover CRISPRs (range: 6-9). Default: 8'),
    OptionParam('-minNR', help='Minimum number of repeats a CRISPR must contain. Default: 3'),
    OptionParam('-minRL', help='Minimum length of the CRISPR repeats. Default: 23'),
    OptionParam('-maxRL', help='Maximum length of the CRISPR repeats. Default: 47'),
    OptionParam('-minSL', help='Minimum length of the CRISPR spacers. Default: 26'),
    OptionParam('-maxSL', help='Maximum length of the CRISPR spacers. Default: 50'),

    OptionParam('-gff',
                help='Output summary results in gff format containing only the positions of the CRISPR arrays. Default: false'),
    OptionParam('-gffFull',
                help='Output detailed results in gff format containing positions of CRISPR arrays and all repeat units. Default: false'),
    OptionParam('-spacers',
                help='Output a fasta formatted file containing the spacers. Default: false'),
    OptionParam('-h', help='Output this handy help message'),
    ArgmntParam(name='query', help='input file of fna sequence'),
    ArgmntParam(name='out', help='output file')]


def run(query, out_dir, gff=True, **kwargs):
    '''Predict CRISPRs for the input file.

    Notes
    -----
    It will create 1 or 2 output files, depending on the parameters:

      1. file containing CRIPSR information, including locations of CRISPRs
         and their sequence composition OR

         GFF file with short information on CRISPR locations OR

         GFFFull file with detailed information on CRISPR locations

      2. (OPTIONAL; -spacers flag) Fasta file of predicted CRISPR spacers

    Parameters
    ----------
    query : str
        input file path of nucleotide sequence
    out_dir : str
        output dir
    gff : bool
        output in gff3 format
    kwargs : dict
        keyword arguments. Other command line parameters for MinCED. key is the option
        (e.g. "-searchWL") and value is the value for the option (e.g. "6").
        If the option is a flag, set the value to None.

    Returns
    -------
    '''
    logger = getLogger(__name__)
    minced = Dumpling('minced', params=Parameters(*params),
                      version='0.2.0', url='https://github.com/ctSkennerton/minced')
    prefix = splitext(basename(query))[0]
    out = join(out_dir, '{}.gff'.format(prefix))
    minced.update(query=query, out=out, gff=gff, **kwargs)
    logger.info('Running {}'.format(minced.command))
    p = minced()
    if not exists(out):
        # minced raises Java error but return 0 exit code if the input file is invalid,
        # so raise the exception for minced manually if no output is created.
        p.returncode = 1
        raise CalledProcessError(
            p.returncode,
            cmd=repr(p.args))
    return minced
