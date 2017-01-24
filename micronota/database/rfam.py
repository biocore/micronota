# ----------------------------------------------------------------------------
# Copyright (c) 2016--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from logging import getLogger

from ..util import split, SplitterTail


logger = getLogger(__name__)


def filter_models(ifile, ofile, negate=False, models={('RF00001', '5S_rRNA'),
                                                      ('RF00002', '5_8S_rRNA'),
                                                      # Permuted mitochondrial genome encoded 5S rRNA
                                                      ('RF02547', 'mtPerm_5S'),
                                                      ('RF01118', 'PK-G12rRNA'),
                                                      ('RF00177', 'SSU_rRNA_bacteria'),
                                                      ('RF01959', 'SSU_rRNA_archaea'),
                                                      ('RF01960', 'SSU_rRNA_eukarya'),
                                                      ('RF02542', 'SSU_rRNA_microsporidia'),
                                                      ('RF02540', 'LSU_rRNA_archaea'),
                                                      ('RF02541', 'LSU_rRNA_bacteria'),
                                                      ('RF02543', 'LSU_rRNA_eukarya'),
                                                      # Trypanosomatid mitochondrial rRNA
                                                      ('RF02545', 'SSU_trypano_mito'),
                                                      ('RF02546', 'LSU_trypano_mito'),
                                                      ('RF00005', 'tRNA'),
                                                      ('RF01852', 'tRNA-Sec'),
                                                      # Mitochondrion encoded tmRNA
                                                      ('RF02544', 'mt_tmRNA'),
                                                      # Alphaproteobacteria transfer messenger RNA
                                                      ('RF01849', 'alpha_tmRNA'),
                                                      # Betaproteobacteria transfer messenger RNA
                                                      ('RF01850', 'beta_tmRNA'),
                                                      # Cyanobacteria transfer messenger RNA
                                                      ('RF01851', 'cyano_tmRNA')}):
    '''Filter away some cm models.

    Parameters
    ----------
    ifile : file-like
        input file of rfam files
    ofile : file-like
        output file with some models filtered away
    negate : bool
        negate the filtering. Keep the specified instead of filtering away them.
    models : Iterable
        list of models to filter away. Default is a list of tRNA, tmRNA, and
        5S/5.8S/16S/18S/23S/28S rRNA
    '''
    splitter = SplitterTail(lambda s: s == '//\n')
    gen = split(splitter)
    j = 0
    i = 0
    for i, record in enumerate(gen(ifile), 1):
        name = record[1].split()[1]
        accn = record[2].split()[1]
        accn_name = (accn, name)
        discard = accn_name in models
        if negate is True:
            discard = not discard
        if discard:
            # logger.debug('Filter %s : %s' % accn_name)
            j += 1
            continue
        else:
            for line in record:
                ofile.write(line)
    logger.debug('Processed %d and filtered %d cm and hmm models' % (i, j))
