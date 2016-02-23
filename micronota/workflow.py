# ----------------------------------------------------------------------------
# Copyright (c) 2015--, micronota development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.workflow import Workflow, method, requires


class AnnotateWF(Workflow):
    '''Annotation workflow.'''
    def initialize_state(self, seq):
        self.state = seq
        self.min_len = 100
        _, self.seq_fp = mkstemp()
        seq.write(self.seq_fp)

    @method(priority=100)
    def check_length(self):
        if len(self.state) < self.min_len:
            self.failed = True

    @method(priority=70)
    @requires(option='cds', values=True)
    def cds(self):
        '''Identify CDS.'''
        print(self.seq_fp)

    @method(priority=60)
    @requires(option='uniref', values=True)
    def search_uniref(self):
        '''Identify homologs and get metadata from them.'''

    @method(priority=50)
    @requires(option='tigrfam', values=True)
    def search_tigrfam(self):
        '''Identify homologs and get metadata from them.'''

    @method(priority=80)
    def ncrna(self):
        '''Identify ncRNA.'''

    @method(priority=90)
    @requires(option='trna', values=True)
    def trna(self):
        '''Identify tRNA and tmRNA.'''


annotate_wf = AnnotateWF(state=None, options={'cds': True})

def annotate_fasta(fp):
    for result in annotate_wf(read(fp, format='fasta', constructor=DNA)):
