__author__ = 'Aurimas Sadauskas'
_version__ = '1.0'
_lastUpdate = '2014-11-23'

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class Analyzer:
    # Constructor
    def __init__(self, program, database, entrez_query, format_type, valid_query_coverage):
        self.program = program
        self.record = ""
        self.database = database
        self.entrez_query = entrez_query
        self.format_type = format_type
        self.valid_query_coverage = valid_query_coverage

    # Method to call blast search, we save result to a file in case we can check what happened
    def startBlastSearch(self, file = "search_output.xml"):
        blast_output = NCBIWWW.qblast(self.program,
            self.database,
            self.record.seq,
            entrez_query=self.entrez_query,
            format_type=self.format_type)
        save_file = open(file, "w")
        save_file.write(blast_output.read())
        return

    # Method to return result of qblast
    # we call quickBlastSearch and specify a file to save xml result to
    def startBlast(self, file = "search_output.xml"):
        self.startBlastSearch(file)

        result_handle = open(file)
        result = NCBIXML.read(result_handle)
        return result;

    # Returns only those records with query coverage over 80 percents
    def blast(self, file = "search_output.xml"):
        brecord = self.startBlast(file)
        not_valid_aligments = []
        query_letters = brecord.query_letters

        # Calculate each alignment query coverage and add invalid records to not_valid_alignments
        for alignment in brecord.alignments:
            match_len_sum = 0
            for hsp in alignment.hsps:
                match_len_sum += hsp.align_length
            val = (match_len_sum*100)/query_letters
            if val < self.valid_query_coverage:
                not_valid_aligments.append(alignment)

        # Remove invalid alignments
        for alignment in not_valid_aligments:
            brecord.alignments.remove(alignment)

        # Return only valid records
        return brecord

    # Method to load main sequence
    def loadMainAlbumin(self, path, format="fasta"):
        self.record = SeqIO.read(open(path), format=format)

    # Method to save blast record in FASTA format
    def saveBlastRecordAsFasta(self, brecord, save_to = "blast_record.fasta"):
        fasta_list = []
        for alignment in brecord.alignments:
            for hsp in alignment.hsps:
                fasta_list.append('>')
                fasta_list.append(alignment.title)
                fasta_list.append('\n')
                fasta_list.append(hsp.sbjct)
                fasta_list.append('\n')
                fasta_list.append('\n')
        with open(save_to,"w") as temp_file:
            temp_file.writelines(fasta_list)

def main():
    part_length_to_analyze = 15

    seqAn = Analyzer("blastp", "swissprot", '"serum albumin"[Protein name] AND mammals[Organism]', "XML",  80);
    seqAn.loadMainAlbumin("serum_albumin_preproprotein.fasta");

    blast = seqAn.blast();
    seqAn.saveBlastRecordAsFasta(blast, "blast_output.fasta");

# Run main
main()