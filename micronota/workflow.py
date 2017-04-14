# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import join, exists, basename, splitext, expanduser
from logging import getLogger
from importlib import import_module
from time import gmtime, strftime

from pkg_resources import resource_filename
from snakemake import snakemake
from skbio import read, write, DNA
import yaml
import numpy as np

from . import module
from .util import _add_cds_metadata
from .quality import compute_gene_score, compute_trna_score, compute_rrna_score, compute_seq_score
from . import __version__


logger = getLogger(__name__)


def annotate(in_fp, in_fmt, min_len, out_dir, out_fmt,
             gcode, kingdom, mode, task,
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
    logger.debug('working dir: %s' % out_dir)
    if force:
        logger.debug('run in force mode, will overwrite existing files')
    if dry_run:
        logger.debug('run in dry mode, will not produce output')

    os.makedirs(out_dir, exist_ok=True)

    out_prefix = join(out_dir, basename(splitext(in_fp)[0]))
    seq_fn_val = out_prefix + '.fna'
    validate_seq(in_fp, in_fmt, min_len, seq_fn_val)

    if config is None:
        config = resource_filename(__package__, kingdom + '.yaml')
    logger.debug('set annotation in %s mode' % mode)
    logger.debug('set annotation as %s' % kingdom)
    logger.debug('use config file: %s' % config)
    with open(config) as fh:
        cfg = yaml.load(fh)

    general = cfg.pop('general', {})
    rules = {}
    if not task:
        task = [i for i in cfg]
    for k, v in cfg.items():
        # specify the annotation task
        if k in task:
            if v is not None:
                for vk, vv in v.items():
                    if vk in rules:
                        raise ValueError('You have multiple config for rule %s' % vk)
                    rules[vk] = vv

    if 'prodigal' in rules:
        param = '%s -g %d' % (rules['prodigal']['params'], gcode)
        if mode == 'finished':
            param = '-p single -c ' + param
        elif mode == 'draft':
            param = '-p single ' + param
        elif mode == 'metagenome':
            param = '-p meta ' + param
        rules['prodigal']['params'] = param
    if 'aragorn' in rules:
        rules['aragorn']['params'] = '%s -gc%d' % (rules['aragorn']['params'], gcode)
    if 'rnammer' in rules:
        rules['rnammer']['params'] = '-S %s %s' % (kingdom[:3], rules['rnammer']['params'])

    logger.debug('run snakemake workflow')
    # only run the targets specified in the yaml file
    targets = list(rules.keys())
    if not targets:
        logger.warning('No annotation task to run')
        return
    rules['seq'] = seq_fn_val

    cfg_file = join(out_dir, 'snakemake.yaml')
    with open(cfg_file, 'w') as out:
        yaml.dump(rules, out, default_flow_style=False)

    snakefile = resource_filename(__package__, 'Snakefile')
    logger.info('Run annotation pipeline')
    success = snakemake(
        snakefile,
        cores=cpus,
        targets=targets,
        # set work dir to output dir so simultaneous runs
        # doesn't interfere with each other.
        workdir=out_prefix,
        printshellcmds=True,
        dryrun=dry_run,
        forceall=force,
        # config=cfg,
        configfile=cfg_file,
        keep_target_files=True,
        # provide this dummy to suppress unnecessary log
        log_handler=lambda s: None,
        quiet=True,  # do not print job info
        keep_logger=False)

    if success:
        # if snakemake finishes successfully
        out_fp = '%s.%s' % (out_prefix, out_fmt)
        protein_xref = general.get('protein_xref')
        if protein_xref is not None:
            protein_xref = expanduser(protein_xref)
        seqs = integrate(seq_fn_val, out_prefix, protein_xref, out_fp, out_fmt)

        # create faa file
        if 'CDS' in task:
            logger.info('Write out protein file')
            faa_fp = out_prefix + '.faa'
            with open(faa_fp, 'w') as out:
                create_faa(seqs.values(), out)

        logger.info('Write summary of the annotation')
        if mode != 'metagenome':
            with open(out_prefix + '.summary', 'w') as out:
                if mode == 'finish':
                    contigs = False
                else:
                    contigs = True
                seq_score = compute_seq_score(seqs.values(), contigs)
                trna_score = rrna_score = gene_score = np.nan
                if 'tRNA' in task:
                    trna_score = compute_trna_score((i.interval_metadata for i in seqs.values()))
                if 'rRNA' in task:
                    rrna_score = compute_rrna_score((i.interval_metadata for i in seqs.values()))
                if 'CDS' in task:
                    gene_score = compute_gene_score(faa_fp)
                out.write('# seq_score: %.2f  tRNA_score: %.2f  rRNA_score: %.2f  gene_score: %.2f\n' % (
                    seq_score, trna_score, rrna_score, gene_score))
                summarize(seqs.values(), out)

    logger.info('Done with annotation')


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
        logger.debug('the sequence file already exists. skip validating step.')
        return
    logger.info('Filter and validate input sequences')
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


def integrate(seq_fp, annot_dir, protein_xref, out_fp, out_fmt='gff3'):
    '''integrate all the annotations and write to disk.

    seq_fn : str
        input seq file name.
    out_dir : str
        annotation output directory.
    out_fmt : str
        output format
    '''
    logger.info('Integrate annotation for output')
    seqs = {}
    for seq in read(seq_fp, format='fasta'):
        seqs[seq.metadata['id']] = seq

    rules = {splitext(f)[0] for f in os.listdir(annot_dir) if f.endswith('.ok')}
    if 'diamond' in rules:
        rules.discard('diamond')
        mod = import_module('.diamond', module.__name__)
        diamond = mod.Module(directory=annot_dir)
        diamond.parse(metadata=protein_xref)
        protein = diamond.result
    else:
        protein = {}
    for rule in rules:
        logger.debug('parse the result from %s output' % rule)
        mod = import_module('.%s' % rule, module.__name__)
        obj = mod.Module(directory=annot_dir)
        obj.parse()
        for seq_id, imd in obj.result.items():
            seq = seqs[seq_id]
            imd._upper_bound = len(seq)
            if rule == 'prodigal':
                cds_metadata = protein.get(seq_id, {})
                _add_cds_metadata(imd, cds_metadata)
            seq.interval_metadata.merge(imd)

    # write out the annotation
    if out_fmt == 'genbank':
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
        write(((sid, seq.interval_metadata) for sid, seq in seqs.items()),
              into=out_fp, format=out_fmt)
    else:
        raise ValueError('Unknown specified output format: %r' % out_fmt)

    return seqs


def summarize(seqs, out):
    '''Summarize the sequences and their annotations.

    Parameters
    ----------
    seqs : list of ``Sequence`` objects
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


def create_faa(seqs, out, genetic_code=11):
    '''Create protein sequence file.

    It creates protein sequences based on the interval features
    with type of "CDS".

    Parameters
    ----------
    seqs : iterable of ``Sequence``
        The list of DNA/RNA sequences
    out : file object
        File object for output
    genetic_code : int
        The fallback genetic code to use
    '''
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

