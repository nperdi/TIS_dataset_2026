import sys
import re
import os

###########################################################
# Input / Output paths
###########################################################

inputDir = "data/ensembl/"

inputFaFileName = os.path.join(
    inputDir,
    "ensembl_proteinCoding_FLANK100_ofProteinCoding_withUNIPROTID.fa"
)

outputDir = "output/positive_protein_coding_representatives/"
os.makedirs(outputDir, exist_ok=True)

###########################################################
# Step 1: FASTA -> TAB
###########################################################

invalidDict = {}

sequenceFaFileName = inputFaFileName
baseName = os.path.splitext(os.path.basename(sequenceFaFileName))[0]

sequenceTabFileName = os.path.join(outputDir, baseName + ".tab")

sequenceFaFile = open(sequenceFaFileName, "r")
sequenceTabFile = open(sequenceTabFileName, "w")

header = ""
for aLine in sequenceFaFile:
    aLine = aLine.strip()
    if len(aLine) == 0:
        continue

    if aLine[0] == ">":
        header = aLine
    else:
        sequenceTabFile.write(header + "\t" + aLine.upper() + "\n")
        if ("N" in aLine) or ("n" in aLine):
            invalidDict[header] = 1

sequenceFaFile.close()
sequenceTabFile.close()

###########################################################
# Step 2: Keep biggest transcript per protein ID
###########################################################

sequenceTabFile = open(sequenceTabFileName, "r")

biggestTranscriptTabFileName = os.path.join(outputDir, baseName + "_biggestTranscript.tab")
biggestTranscriptFaFileName = os.path.join(outputDir, baseName + "_biggestTranscript.fa")
biggestTranscriptBedFileName = os.path.join(outputDir, baseName + "_biggestTranscript.bed")
biggestTranscriptWithoutFlankBedFileName = os.path.join(
    outputDir, baseName + "_biggestTranscript_withoutFlank.bed"
)

biggestTranscriptTabFile = open(biggestTranscriptTabFileName, "w")
biggestTranscriptFaFile = open(biggestTranscriptFaFileName, "w")
biggestTranscriptBedFile = open(biggestTranscriptBedFileName, "w")
biggestTranscriptWithoutFlankBedFile = open(biggestTranscriptWithoutFlankBedFileName, "w")

proteinDict = {}

for aLine in sequenceTabFile:
    # Example header:
    # >ENSG00000154227|15|Q8IU89|-1|100400395|100544720|ENST00000284382
    aLine = aLine.rstrip().strip(" ")
    if len(aLine) == 0:
        continue

    cols = aLine.split("\t")
    if len(cols) < 2:
        continue

    header = cols[0]
    sequence = cols[1]

    headerList = header[1:].split("|")
    if len(headerList) < 7:
        continue

    uniprotId = headerList[2]
    transcriptLength = len(sequence)

    match = re.search('[ACTGactg]', sequence)
    if match:
        if uniprotId in proteinDict:
            storedTranscriptLen = proteinDict[uniprotId][1]
            if int(storedTranscriptLen) < transcriptLength:
                proteinDict[uniprotId] = [aLine, transcriptLength]
        else:
            proteinDict[uniprotId] = [aLine, transcriptLength]

sequenceTabFile.close()

print("len(proteinDict)", len(proteinDict))

###########################################################
# Step 3: Write representative transcript outputs
###########################################################

for aUniprotId in proteinDict:
    cols = proteinDict[aUniprotId][0].split("\t")
    header = cols[0]
    sequence = cols[1]

    faRecord = header + "\n" + sequence + "\n"

    headerFields = header.split("|")
    strandNum = headerFields[3]

    chrom = "chr" + headerFields[1]
    start = headerFields[4]
    end = headerFields[5]
    transcriptId = headerFields[6]
    geneId = headerFields[0][1:]
    proteinId = headerFields[2]

    if strandNum == "1":
        strand = "+"
        removeFlankingStart = str(int(start) + 100)
        removeFlankingEnd = end
    else:
        strand = "-"
        removeFlankingStart = start
        removeFlankingEnd = str(int(end) - 100)

    bedName = geneId + "|" + transcriptId + "|" + proteinId

    fullBedLine = "\t".join([
        chrom,
        start,
        end,
        bedName,
        "0",
        strand
    ]) + "\n"

    withoutFlankBedLine = "\t".join([
        chrom,
        removeFlankingStart,
        removeFlankingEnd,
        bedName,
        "0",
        strand
    ]) + "\n"

    biggestTranscriptTabFile.write(proteinDict[aUniprotId][0] + "\n")
    biggestTranscriptFaFile.write(faRecord)
    biggestTranscriptBedFile.write(fullBedLine)
    biggestTranscriptWithoutFlankBedFile.write(withoutFlankBedLine)

biggestTranscriptTabFile.close()
biggestTranscriptFaFile.close()
biggestTranscriptBedFile.close()
biggestTranscriptWithoutFlankBedFile.close()

print("Output directory:", outputDir)
print("TAB :", biggestTranscriptTabFileName)
print("FA  :", biggestTranscriptFaFileName)
print("BED :", biggestTranscriptBedFileName)
print("BED without flank:", biggestTranscriptWithoutFlankBedFileName)