# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, skipIf, main
from pkg_resources import resource_filename
from shutil import which
from os.path import join
from os import mkdir
from tempfile import TemporaryDirectory, mkdtemp

from snakemake.logging import logger, setup_logger
from snakemake import snakemake


class TestRules(TestCase):
    def setUp(self):
        self.seq_fn = 'test.fna'
        self.snakefile = resource_filename('micronota', 'rules/Snakefile')
        self.config = {'seq': 'test.fna', 'genetic_code': 11}
        self.tmpd = mkdtemp()

    def _run_snakemake(self, config):
        open(join(self.tmpd, self.seq_fn), 'w').close()
        success = snakemake(
            self.snakefile,
            cores=8,
            workdir=self.tmpd,
            dryrun=True,
            forcetargets=True,
            config=self.config,
            keep_target_files=True,
            keep_logger=True)
        with open(logger.logfile) as o:
            return o.read()

    @skipIf(which("aragorn") is None, 'Aragorn not installed.')
    def test_aragorn(self):
        self.config['tools'] = {'aragorn':
                                {'params': '-l',
                                 'priority': 50,
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('aragorn -l -gc11 -o aragorn/{0}.txt {0} &> aragorn/{0}.log && '
               'touch aragorn/{0}.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("minced") is None, 'minced not installed.')
    def test_minced(self):
        self.config['tools'] = {'minced':
                                {'params': '',
                                 'priority': 50,
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('minced  -gff {0} minced/{0}.gff &> minced/{0}.log && '
               'touch minced/{0}.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("prodigal") is None, 'Prodigal not installed.')
    def test_prodigal(self):
        self.config['tools'] = {'prodigal':
                                {'params': '-p meta',
                                 'priority': 90,
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('prodigal -p meta -f gff -g 11 -i {0} -o prodigal/{0}.gff'
               ' -a prodigal/{0}.faa -d prodigal/{0}.fna &>'
               ' prodigal/{0}.log && '
               'touch prodigal/{0}.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("cmscan") is None, 'Infernal not installed.')
    def test_cmscan(self):
        with TemporaryDirectory() as tmpd:
            db = join(tmpd, 'rfam.cm')
            open(db, 'w').close()
            self.config['tools'] = {'cmscan':
                                    {'params': '',
                                     'priority': 30,
                                     'threads': 2,
                                     'db': db}}
            log = self._run_snakemake(self.config)
        exp = ('cmscan  --cpu {n} --tblout cmscan/{o}.txt {db} {o} &>'
               ' cmscan/{o}.log && touch cmscan/{o}.ok').format(
                   o=self.seq_fn, db=db, n=2)
        self.assertIn(exp, log)

    @skipIf(which("diamond") is None, 'diamond not installed.')
    def test_diamond(self):
        with TemporaryDirectory() as tmpd:
            db = tmpd
            open(join(db, 'uniref90.dmnd'), 'w').close()
            open(join(db, 'uniref50.dmnd'), 'w').close()
            d = join(self.tmpd, 'prodigal')
            mkdir(d)
            open(join(d, self.seq_fn + '.faa'), 'w').close()
            config = {'params': '--index-chunks 1 --query-cover 80 -k 3',
                      'priority': 80, 'threads': 4, 'db': db}
            self.config['tools'] = {'diamond': config}
            log = self._run_snakemake(self.config)

        exp = ('diamond blastp {p} --threads {n} --db {db}/uniref90.dmnd -q {i}'
               ' -a {o}_uniref90.daa &> {o}_uniref90.daa.log'.format(
                   p=config['params'],
                   n=config['threads'],
                   db=config['db'],
                   o='diamond/%s' % self.seq_fn,
                   i='prodigal/%s.faa' % self.seq_fn))
        self.assertIn(exp, log)

        unmatched = join('diamond', self.seq_fn + '_uniref90_unmatched.faa')
        exp = ('diamond blastp {p} --threads {n} --db {db}/uniref50.dmnd -q {i}'
               ' -a {o}_uniref50.daa &> {o}_uniref50.daa.log'.format(
                   p=config['params'],
                   n=config['threads'],
                   db=config['db'],
                   o='diamond/%s' % self.seq_fn, i=unmatched))
        self.assertIn(exp, log)


if __name__ == '__main__':
    # this is required to set snakemake logging to files correctly
    setup_logger(printshellcmds=True)
    logger.logger.removeHandler(logger.logger.handlers[1])

    main()
