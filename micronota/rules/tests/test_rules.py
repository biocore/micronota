# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, skipIf, main
from shutil import which, rmtree
from os.path import join, abspath, dirname
from os import mkdir
from tempfile import mkdtemp
from subprocess import Popen, PIPE


class TestRules(TestCase):
    def setUp(self):
        self.tmpd = mkdtemp()
        self.snakefile = join(self.tmpd, 'Snakefile')
        self.path = abspath(dirname(__file__))
        self.seq_fn = 'test.fna'
        self.out_path = join(self.tmpd, self.seq_fn)
        self.template = '''workdir: '{wd}'
fn = 'test.fna'
tool_output = '{o}'
tool_config = %r
config['seq'] = '{i}'
config['genetic_code'] = 2
include: '%s'
'''.format(wd=self.tmpd, i=self.input_seq, o=self.out_path)

    def _run_snakemake(self, tool, s):
        with open(self.snakefile, 'w') as out:
            out.write(self.template % (
                s, join(dirname(self.path), '%s.rule' % tool)))
        p = Popen(['snakemake', '-s', self.snakefile, '-F', '-p', '-n',
                   '--cores', '8',
                   '--keep-target-files', self.out_path + '.gff',
                   '--config', 'output_dir=%s' % self.tmpd],
                  universal_newlines=True, stdout=PIPE, )
        stdout, stderr = p.communicate()
        return stdout

    @skipIf(which("aragorn") is None, 'Aragorn not installed.')
    def test_aragorn(self):
        stdout = self._run_snakemake(
            'aragorn', {'params': '-l', 'priority': 50, 'threads': 1})
        exp = 'aragorn -l -gc2 -o {o}.txt {i} &> {o}.log'.format(
            o=self.out_path, i=self.input_seq)
        self.assertIn(exp, stdout)
        self.assertIn('priority: 50', stdout)

    @skipIf(which("minced") is None, 'minced not installed.')
    def test_minced(self):
        stdout = self._run_snakemake(
            'minced', {'params': '', 'priority': 50, 'threads': 1})
        exp = 'minced  -gff {i} {o}.gff &> {o}.log'.format(
            o=self.out_path, i=self.input_seq)
        self.assertIn(exp, stdout)
        self.assertIn('priority: 50', stdout)

    @skipIf(which("prodigal") is None, 'Prodigal not installed.')
    def test_prodigal(self):
        stdout = self._run_snakemake(
            'prodigal', {'params': '-p meta', 'priority': 90, 'threads': 1})
        exp = ('prodigal -p meta -f gff -g 2 -i {i} -o {o}.gff'
               ' -a {o}.faa -d {o}.fna &> {o}.log'.format(
                   o=self.out_path, i=self.input_seq))
        self.assertIn(exp, stdout)
        self.assertIn('priority: 90', stdout)

    @skipIf(which("cmscan") is None, 'Infernal not installed.')
    def test_cmscan(self):
        db = join(self.tmpd, 'rfam.cm')
        open(db, 'w').close()
        config = {'params': '', 'priority': 30, 'threads': 2, 'db': db}
        stdout = self._run_snakemake('cmscan', config)
        exp = 'cmscan  --cpu {n} --tblout {o}.txt {db} {i} &> {o}.log'.format(
            o=self.out_path, i=self.input_seq, db=config['db'], n=config['threads'])
        self.assertIn(exp, stdout)
        self.assertIn('priority: 30', stdout)

    @skipIf(which("diamond") is None, 'diamond not installed.')
    def test_diamond(self):
        db = self.tmpd
        open(join(db, 'uniref90.dmnd'), 'w').close()
        open(join(db, 'uniref50.dmnd'), 'w').close()
        mkdir(join(self.tmpd, 'prodigal'))
        i = join(self.tmpd, 'prodigal', self.seq_fn + '.faa')
        open(i, 'w').close()
        config = {'params': '--index-chunks 1 --query-cover 80 -k 3',
                  'priority': 80, 'threads': 4, 'db': db}
        stdout = self._run_snakemake('diamond', config)
        exp = ('diamond blastp {p} --threads {n} --db {db}/uniref90.dmnd -q {i}'
               ' -a {o}_uniref90.daa &> {o}_uniref90.daa.log'.format(
                   p=config['params'],
                   n=config['threads'],
                   db=config['db'],
                   o=self.out_path, i=i))
        self.assertIn(exp, stdout)
        unmatched = join(self.tmpd, self.seq_fn + '_uniref90_unmatched.faa')
        exp = ('diamond blastp {p} --threads {n} --db {db}/uniref50.dmnd -q {i}'
               ' -a {o}_uniref50.daa &> {o}_uniref50.daa.log'.format(
                   p=config['params'],
                   n=config['threads'],
                   db=config['db'],
                   o=self.out_path, i=unmatched))
        self.assertIn(exp, stdout)
        self.assertIn('priority: 80', stdout)

    def tearDown(self):
        rmtree(self.tmpd)


if __name__ == '__main__':
    main()
