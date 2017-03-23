# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from collections import defaultdict
from os.path import join, exists, basename, splitext
from logging import getLogger
from importlib import import_module
from time import gmtime, strftime

from pkg_resources import resource_filename
from snakemake import snakemake
from skbio import read, write, DNA
import yaml

from . import parsers
from .parsers.cds import _add_cds_metadata, _fetch_cds_metadata
from .quality import compute_gene_score, compute_trna_score, compute_rrna_score, compute_seq_score
from . import __version__


logger = getLogger(__name__)


def annotate(in_fp, in_fmt, min_len, out_dir, out_fmt, gcode,
             mode, kingdom,
             cpus, force, dry_run, config):
    '''Annotate the sequences in the input file.

    Parameters
    ----------
    in_fp : str
        Input seq file name
    in_fmt : str
        Input file format
    min_len : int
        The threshold of seq length to be filtered away
    out_dir : str
        Output file directory.
    out_fmt : str
        Output file format
    gcode : int
        The translation table to use for protein-coding genes
    mode : bool
        Run with metagenomic mode?
    kingdom : str
        The kingdom where the sequences are from
    cpus : int
        Number of cpus to use.
    force : bool
        Force to overwrite.
    dry_run : bool
    config : config file for snakemake
    '''
    logger.info('Running annotation pipeline')

    os.makedirs(out_dir, exist_ok=True)

    seq_fn = basename(in_fp)

    seq_fn_val = splitext(seq_fn)[0] + '.valid.fna'
    validate_seq(in_fp, in_fmt, min_len, join(out_dir, seq_fn_val))

    if config is None:
        config = resource_filename(__package__, 'config.yaml')
    with open(config) as fh:
        cfg = yaml.load(fh)

    cfg['seq'] = seq_fn_val
    cfg['genetic_code'] = gcode
    cfg['kingdom'] = kingdom
    prodigal = cfg['structural_annotation']['prodigal']
    if mode == 'finished':
        prodigal['params'] = '-p single -c'
    elif mode == 'draft':
        prodigal['params'] = '-p single'
    elif mode == 'metagenome':
        prodigal['params'] = '-p meta'
    cfg['mode'] = mode
    cfg_file = join(out_dir, 'config.yaml')
    with open(cfg_file, 'w') as out:
        yaml.dump(cfg, out, default_flow_style=False)

    snakefile = resource_filename(__package__, 'rules/Snakefile')

    success = snakemake(
        snakefile,
        cores=cpus,
        # set work dir to output dir so simultaneous runs
        # doesn't interfere with each other.
        workdir=out_dir,
        printshellcmds=True,
        dryrun=dry_run,
        forceall=force,
        # config=cfg,
        configfile=cfg_file,
        keep_target_files=True,
        keep_logger=False)

    if success:
        # if snakemake finishes successfully
        integrate(cfg, out_dir, seq_fn_val, out_fmt)

    logger.info('Done with annnotation')


def validate_seq(in_fp, in_fmt, min_len, out_fp):
    '''Validate input seq file.

    1. filter out short seq;
    2. validate seq IDs (no duplicates)
    3. remove gaps in the sequence if there is any
    4. convert to fasta format

    Parameters
    ----------
    in_fp : str
        input seq file path
    in_fmt : str
        the format of seq file
    min_len : int
        cutoff of seq len to filter away
    out_fp : str
        output seq file path
    '''
    if exists(out_fp):
        # do not overwrite because all the snakemake steps will be rerun when
        # this file is updated.
        logger.debug('The sequence file already exist. Skip validating step.')
        return
    logger.info('Filtering and validating input sequences')
    ids = set()
    with open(out_fp, 'w') as out:
        for seq in read(in_fp, format=in_fmt, constructor=DNA):
            seq = seq.degap()
            if len(seq) < min_len:
                continue
            if in_fmt == 'genbank':
                seq.metadata['id'] = seq.metadata['LOCUS']['locus_name']
            try:
                ident = seq.metadata['id']
            except KeyError:
                raise KeyError('Ill input file format: at least one sequences do not have IDs.')
            if ident in ids:
                raise ValueError(
                    'Duplicate seq IDs in your input file: {}'.format(ident))
            else:
                ids.add(ident)
            write(seq, format='fasta', into=out)


def integrate(cfg, out_dir, seq_fn, out_fmt='genbank'):
    '''integrate all the annotations and write to disk.

    seq_fn : str
        input seq file name.
    out_dir : str
        annotation output directory.
    out_fmt : str
        output format
    '''
    logger.info('Integrating annotation for output')
    seqs = {}
    for seq in read(join(out_dir, seq_fn), format='fasta'):
        seqs[seq.metadata['id']] = seq

    hits = []
    for f in os.listdir(out_dir):
        if f.endswith('.hit'):
            # protein homologous map
            hits.append(join(out_dir, f))
        elif f.endswith('.ok'):
            try:
                # parse the output of each feature prediction tool into
                # interval metadata
                tool = f.split('.', 1)[0]
                submodule = import_module('.%s' % tool, parsers.__name__)
                func = getattr(submodule, 'parse')
                fp = join(out_dir, f.replace('.ok', '.txt'))
                for seq_id, imd in func(fp):
                    seq = seqs[seq_id]
                    imd._upper_bound = len(seq)
                    seq.interval_metadata.merge(imd)
            except Exception as e:
                # print more useful info if exception is raised
                raise Exception('Error while parsing %s into ``IntervalMetadata``' % f) from e
    # add functional metadata to the protein-coding gene
    # create defaultdict of defaultdict of dict
    cds_metadata = defaultdict(lambda: defaultdict(dict))
    for f in hits:
        for seq_id, idx, md in _fetch_cds_metadata(f, cfg['metadata']):
            cds_metadata[seq_id][idx].update(md)
    if cds_metadata:
        _add_cds_metadata(seqs, cds_metadata)

    # write out the annotation
    if out_fmt == 'genbank':
        out_fp = join(out_dir, '%s.gbk' % seq_fn)
        with open(out_fp, 'w') as out:
            for sid, seq in seqs.items():
                seq.metadata['LOCUS'] = {
                    'locus_name': sid,
                    'size': len(seq),
                    'unit': 'bp',
                    'mol_type': 'DNA',
                    'shape': 'linear',
                    'division': None,
                    'date': strftime("%d-%b-%Y", gmtime())}
                seq.metadata['ACCESSION'] = ''
                seq.metadata['VERSION'] = ''
                seq.metadata['KEYWORDS'] = '.'
                seq.metadata['SOURCE'] = {'ORGANISM': 'genus species', 'taxonomy': 'unknown'}
                seq.metadata['COMMENT'] = 'Annotated with %s %s' % (__package__, __version__)
                write(seq, into=out, format=out_fmt)
    elif out_fmt == 'gff3':
        out_fp = join(out_dir, '%s.gff' % seq_fn)
        write(((sid, seq.interval_metadata) for sid, seq in seqs.items()),
              into=out_fp, format=out_fmt)
    else:
        raise ValueError('Unknown specified output format: %r' % out_fmt)

    # create faa file
    faa_fp = join(out_dir, '%s.faa' % seq_fn)
    create_faa(seqs.values(), faa_fp)
    if cfg['mode'] != 'metagenome':
        with open(join(out_dir, 'summary.txt'), 'w') as out:
            if cfg['mode'] == 'finish':
                contigs = False
            else:
                contigs = True
            seq_score = compute_seq_score(seqs.values(), contigs)
            trna_score = compute_trna_score((i.interval_metadata for i in seqs.values()))
            rrna_score = compute_rrna_score((i.interval_metadata for i in seqs.values()))
            gene_score = compute_gene_score(faa_fp)
            out.write('# seq_score: %.2f  tRNA_score: %.2f  rRNA_score: %.2f  gene_score: %.2f\n' % (
                seq_score, trna_score, rrna_score, gene_score))
            summarize(seqs.values(), out)


def summarize(seqs, out):
    '''Summarize the sequences and their annotation.

    Parameters
    ----------
    seqs : list of ``Sequence``
    out : file object
        the file object to output to
    '''
    types = ['CDS', 'ncRNA', 'rRNA', 'tRNA',
             'tandem_repeat', 'terminator', 'CRISPR']

    out.write('#seq_id\tlength\tnuc_freq\t')
    out.write('\t'.join(types))
    out.write('\n')
    for seq in seqs:
        freq = seq.frequencies(relative=True)
        items = [seq.metadata['id'], str(len(seq)),
                 ';'.join(['%s:%.2f' % (k, freq[k]) for k in sorted(freq)])]
        imd = seq.interval_metadata
        for t in types:
            feature = imd.query(metadata={'type': t})
            items.append(str(len([i for i in feature])))
        out.write('\t'.join(items))
        out.write('\n')


def create_faa(seqs, out_fp, genetic_code=11):
    '''Create protein sequence file.

    It creates protein sequences based on the interval features
    with type of "CDS".

    Parameters
    ----------
    seqs : iterable of ``Sequence``
        The list of DNA/RNA sequences
    out_fp : str
        File path for output
    genetic_code : int
        The fallback genetic code to use
    '''
    with open(out_fp, 'w') as out:
        for seq in seqs:
            for cds in seq.interval_metadata.query(metadata={'type': 'CDS'}):
                fna = DNA.concat([seq[start:end] for start, end in cds.bounds])
                if cds.metadata.get('strand', '.') == '-':
                    fna = fna.reverse_complement()
                try:
                    # if translation table is not available in metadata, fallback
                    # to what is specified in the func parameter
                    faa = fna.translate(cds.metadata.get('transl_table', genetic_code))
                    faa.metadata['description'] = cds.metadata.get('product', '')
                    # CDS metadata must have key of 'ID'
                    faa.metadata['id'] = cds.metadata['ID']
                    write(faa, into=out, format='fasta')
                except NotImplementedError:
                    logger.warning('This gene has degenerate nucleotide and will not be translated.')
