# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, skipIf, main
from shutil import which, rmtree
from os.path import join
from tempfile import TemporaryDirectory, mkdtemp

from pkg_resources import resource_filename
from snakemake.logging import logger, setup_logger
from snakemake import snakemake


class TestRules(TestCase):
    def setUp(self):
        self.seq_fn = 'test.fna'
        self.snakefile = resource_filename('micronota', 'rules/Snakefile')
        self.tmpd = mkdtemp()
        self.config = {'general': {'seq': 'test.fna', 'genetic_code': 11, 'kingdom': 'bacteria'},
                       'structural_annotation': {},
                       'bacteria': {},
                       'protein': {}}
        # this is required to set snakemake logging to files correctly
        setup_logger(printshellcmds=True)
        logger.logger.removeHandler(logger.logger.handlers[1])

    def _run_snakemake(self, config):
        open(join(self.tmpd, self.seq_fn), 'w').close()
        snakemake(
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
        self.config['structural_annotation'] = {'aragorn':
                                {'params': '-l',
                                 'priority': 50,
                                 'output': 'aragorn',
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('aragorn -l -gc11 -o aragorn.txt {0} &> aragorn.log && '
               'touch aragorn.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("minced") is None, 'minced not installed.')
    def test_minced(self):
        self.config['structural_annotation'] = {'minced':
                                {'params': '',
                                 'output': 'minced',
                                 'priority': 50,
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('minced  -gff {0} minced.gff &> minced.log && '
               'touch minced.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("prodigal") is None, 'Prodigal not installed.')
    def test_prodigal(self):
        self.config['structural_annotation'] = {'prodigal':
                                {'params': '-p meta',
                                 'output': 'prodigal',
                                 'priority': 90,
                                 'threads': 1}}
        log = self._run_snakemake(self.config)
        exp = ('prodigal -p meta -f gff -g 11 -i {0} -o prodigal.txt'
               ' -a prodigal.faa -d prodigal.fna &>'
               ' prodigal.log && '
               'touch prodigal.ok').format(self.seq_fn)
        self.assertIn(exp, log)

    @skipIf(which("cmscan") is None, 'Infernal not installed.')
    def test_cmscan(self):
        with TemporaryDirectory() as tmpd:
            db = join(tmpd, 'rfam.cm')
            open(db, 'w').close()
            self.config['structural_annotation'] = {'cmscan':
                                    {'params': '',
                                     'output': 'cmscan',
                                     'priority': 30,
                                     'threads': 2,
                                     'db': db}}
            log = self._run_snakemake(self.config)
        exp = ('cmscan  --cpu {n} --tblout cmscan.txt {db} {o} &>'
               ' cmscan.log && touch cmscan.ok').format(
                   o=self.seq_fn, db=db, n=2)
        self.assertIn(exp, log)

    def tearDown(self):
        logger.cleanup()
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
