__author__ = 'Aurimas Sadauskas'
_version__ = '1.0'
_lastUpdate = '2014-11-23'

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class Analyzer:
    # Constructor
    def __init__(self, program):
        self.program = program


def main():
    program = "blastp"
    sequence_analyzer = Analyzer(program)
    sequence_analyzer.programName();
# Run main
main()