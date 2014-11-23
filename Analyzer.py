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
    def startBlastSearch(self, input="search_output.xml"):
        blast_output = NCBIWWW.qblast(self.program, self.database, self.record.seq,
            entrez_query = self.entrez_query, format_type = self.format_type)

        output_file = open(input, "w")
        output_file.write(blast_output.read())
        return

    # Method to return result of qblast
    # we call quickBlastSearch and specify a file to save xml result to
    def startBlast(self, input="search_output.xml"):
        self.startBlastSearch(input)

        result = open(input)
        return NCBIXML.read(result);

    # Returns only those records with query coverage over 80 percents
    def blast(self, input="search_output.xml"):
        blastrecords = self.startBlast(input)
        not_valid_aligments = []
        query_letters = blastrecords.query_letters

        # Calculate each alignment query coverage and add invalid records to not_valid_alignments
        for alignment in blastrecords.alignments:
            match_len_sum = 0
            for hsp in alignment.hsps:
                match_len_sum += hsp.align_length
            val = (match_len_sum*100)/query_letters
            if val < self.valid_query_coverage:
                not_valid_aligments.append(alignment)

        # Remove invalid alignments
        for alignment in not_valid_aligments:
            blastrecords.alignments.remove(alignment)

        # Return only valid records
        return blastrecords

    # Method to load main sequence
    def loadMainAlbumin(self, path, format="fasta"):
        self.record = SeqIO.read(open(path), format=format)

    # Method to save blast record in FASTA format
    def saveBlastRecordAsFasta(self, blastrecords, saveFile="blast_record.fasta"):
        fastalist = []
        for alignment in blastrecords.alignments:
            for hsp in alignment.hsps:
                fastalist.append('>')
                fastalist.append(alignment.title)
                fastalist.append('\n')
                fastalist.append(hsp.sbjct)
                fastalist.append('\n')
                fastalist.append('\n')
        with open(saveFile, "w") as tempfile:
            tempfile.writelines(fastalist)

    def startMafft(self, input="blast_record.fasta", output="mafft_output.fasta"):
        os.system('mafft --quiet ' + input + ' > ' + output)
        return

    # Method to calculate mistakes in a part
    def mistakes_in_part(self, main_sequence, all_seqs, start_index, part_length):
        mistakes_count = 0;

        # have to skip the first sequence because we compare with it
        # for all sequences except first
        for i in range(1, len(all_seqs)):
            # take a sequence
            sequence = all_seqs[i]
            # count mistakes in main_sequence compared to that sequence that we took in a range
            for j in range (start_index, start_index + part_length):
                if (sequence[j] != main_sequence[j]):
                    mistakes_count += 1
        return mistakes_count

    # Method to analyze sequence and find the most similar and dissimilar parts
    # main_sequence to compare with must be first in a source file
    def analysis(self, part_length, with_file = "mafft.fasta"):
        # Reading source file to get sequences
        all_seqs = []
        for r in SeqIO.parse(open(with_file), "fasta"):
            all_seqs.append(r.seq)
            # NOTE: `source` file's first line should be the sequence that will be compared with others
        # Let's mark main sequence
        main_sequence = all_seqs[0]

        # initialise markers for the most similar and dissimilar parts
        min_mistakes = -1
        min_index = 0
        max_mistakes = -1
        max_index  = 0

        for i in range(0, len(main_sequence) - part_length):
            # count mistakes starting from i to i+part_length
            mistakes = self.mistakes_in_part(main_sequence, all_seqs, i, part_length)

            # mark max and min mistakes in a part
            if (mistakes > max_mistakes) or (min_mistakes == -1):
                min_mistakes = mistakes
                min_index = i
            if (mistakes < min_mistakes) or (max_mistakes == -1):
                max_mistakes = mistakes
                max_index = i


        print("Maziausiai panasi seka: " +\
              str(min_index) + ".." + str(min_index + part_length) +\
              " (klaidu: " + str(min_mistakes) + ")")

        # Sequence checking
        for sequence in all_seqs:
            print(sequence[min_index:(min_index+part_length)])
            break

        print("Labiausiai panasi seka: " +\
              str(max_index) + ".." + str(max_index + part_length) +\
              " (klaidu: " + str(max_mistakes) + ")")
        # Sequence checking
        for sequence in all_seqs:
            print(sequence[max_index:(max_index+part_length)])
            break
        return

def main():
    seqAn = Analyzer("blastp", "swissprot", '"serum albumin"[Protein name] AND mammals[Organism]', "XML",  80);
    seqAn.loadMainAlbumin("serum_albumin_preproprotein.fasta");

    #blast = seqAn.blast();
    #seqAn.saveBlastRecordAsFasta(blast, "blast_output.fasta");
    #seqAn.startMafft("blast_output.fasta");
    seqAn.analysis(15);

# Run main
main()