__author__ = 'Aurimas Sadauskas'
_version__ = '1.0'
_lastUpdate = '2014-11-24'

import os
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class Analyzer:
    # Constructor
    def __init__(self, program, database, entrezQuery, formatType, queryCoverage):
        self.program = program
        self.record = ""
        self.database = database
        self.entrezQuery = entrezQuery
        self.formatType = formatType
        self.queryCoverage = queryCoverage

    def startBlastSearch(self, input="search_output.xml"):
        blastOutput = NCBIWWW.qblast(self.program, self.database, self.record.seq,
            entrez_query = self.entrezQuery, format_type = self.formatType)

        outputFile = open(input, "w")
        outputFile.write(blastOutput.read())
        return;

    def startBlast(self, input="search_output.xml"):
        self.startBlastSearch(input)

        result = open(input)
        return NCBIXML.read(result);

    def blast(self, input="search_output.xml"):
        blastRecords = self.startBlast(input)
        notValidAlignments = []
        queryLetters = blastRecords.query_letters

        for alignment in blastRecords.alignments:
            sum = 0
            for hsp in alignment.hsps:
                sum += hsp.align_length
            val = (sum * 100) / queryLetters
            if val < self.queryCoverage:
                notValidAlignments.append(alignment)

        for alignment in notValidAlignments:
            blastRecords.alignments.remove(alignment)

        return blastRecords;

    def loadMainAlbumin(self, path, format="fasta"):
        self.record = SeqIO.read(open(path), format=format)

    def saveBlastRecord(self, blastRecords, saveFile="blast_output.fasta"):
        fastaList = []
        for alignment in blastRecords.alignments:
            for hsp in alignment.hsps:
                fastaList.append(">{0}\n{1}\n\n".format(alignment.title, hsp.sbjct))

        with open(saveFile, "w") as tempFile:
            tempFile.writelines(fastaList)

    def startMafft(self, input="blast_record.fasta", output="mafft_output.fasta"):
        os.system("mafft --quiet {0} > {1}".format(input, output))
        return;

    def mistakesInPart(self, mainSequence, allSeqs, startIndex, partLength):
        mistakes = 0

        for i in range(1, len(allSeqs)):
            sequence = allSeqs[i]

            for j in range (startIndex, startIndex + partLength):
                if (sequence[j] != mainSequence[j]):
                    mistakes += 1

        return mistakes;

    def analysis(self, partLength, input="mafft_output.fasta", format="fasta"):
        allSeqs = []
        for r in SeqIO.parse(open(input), format=format):
            allSeqs.append(r.seq)

        mainSequence = allSeqs[0]
        minMistakes = -1
        minIndex = 0
        maxMistakes = -1
        maxIndex  = 0

        for i in range(0, len(mainSequence) - partLength):
            mistakes = self.mistakesInPart(mainSequence, allSeqs, i, partLength)

            if (mistakes > maxMistakes) or (minMistakes == -1):
                minMistakes = mistakes
                minIndex = i
            if (mistakes < minMistakes) or (maxMistakes == -1):
                maxMistakes = mistakes
                maxIndex = i

        print("Mažiausiai panaši seka: {0}".format(allSeqs[0][minIndex:(minIndex+partLength)]))
        print("Labiausiai panaši seka: {0}".format(allSeqs[0][maxIndex:(maxIndex+partLength)]))
        return;

def main():
    searchOutputFile = "search_output.xml"
    blastOutputFile = "blast_output.fasta"

    seqAn = Analyzer("blastp", "swissprot", '"serum albumin"[Protein name] AND mammals[Organism]', "XML",  80)
    seqAn.loadMainAlbumin("serum_albumin_preproprotein.fasta")

    if (os.path.isfile(searchOutputFile) == False):
        blast = seqAn.blast(searchOutputFile)
        seqAn.saveBlastRecord(blast, blastOutputFile)
    if (os.path.isfile(blastOutputFile) == True):
        seqAn.startMafft(blastOutputFile)

    seqAn.analysis(15)

# Run main
main()