"""
Script purpose:
This script extracts candidate intergenic ORFs from previously identified
intergenic genomic regions.

Main steps:
1. Read a BED file containing filtered intergenic regions.
2. Keep up to a fixed number of samples from these regions.
3. For each region, extract a genomic subregion from the reference genome
   (from relativeStart to relativeEnd relative to the region start).
4. Save these extracted regions in BED format.
5. Use bedtools getfasta to retrieve the corresponding DNA sequences
   from the reference genome.
6. Scan each extracted sequence for ATG start codons.
7. For each ATG, search for the first in-frame stop codon (TAA, TAG, TGA).
8. Keep only ORF-like sequences that:
   - have a valid in-frame stop codon,
   - contain enough upstream sequence for the selected flank size,
   - do not contain ambiguous base 'N'.
9. Repeat ORF extraction for different upstream flank sizes
   (for example 0 nt and 100 nt before the ATG).
10. Save the resulting candidate ORFs in:
   - FASTA format
   - BED format

Inputs:
- BED file with filtered intergenic regions
- Reference genome FASTA file

Outputs:
- BED file with extracted intergenic regions
- FASTA file with extracted intergenic region sequences
- FASTA files with candidate intergenic ORFs for each flank size
- BED files with candidate intergenic ORFs for each flank size

Notes:
- The script searches only on the '+' strand as defined in the generated BED entries.
- ORFs are defined here as sequences starting with ATG and ending at the first
  downstream in-frame stop codon.
- The number of extracted regions and ORFs is limited by numOfsamples.
"""


import sys
sys.path.append('../lib')
# from nMersDB import *
import os

genome_file = "../data/GENOME/hg38.fa"
annotationFileName = "output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir = "output/intergenic_ORFs/"
command = "mkdir -p " + workingDir
os.system(command)

###########################################################


###########################################################
numOfsamples = 15000
relativeStart = 0
relativeEnd = 30000
flankRegionList = [0, 100]
stopCodons = {"TAA", "TAG", "TGA"}
###########################################################

annotationFile = open(annotationFileName, "r")

intergenicFlankBedFileName = workingDir + "intergenicRegions.bed"
intergenicFlankBedFile = open(intergenicFlankBedFileName, "w")

cnt = 0
try:
    for aLine in annotationFile:
        if len(aLine) <= 1:
            continue

        cols = aLine.rstrip().split("\t")
        chrom = "chr" + cols[0][3:]
        start = int(cols[1])
        end = int(cols[2])
        name = cols[3]

        currentRelativeEnd = relativeEnd
        if end - start < currentRelativeEnd:
            currentRelativeEnd = end - start

        extractRegionStart = start + relativeStart
        extractRegionEnd = start + currentRelativeEnd

        extractRegionExtendedBedFileLine = "\t".join(
            map(
                str,
                [
                    chrom,
                    extractRegionStart,
                    extractRegionEnd,
                    name + ":" + chrom + ":" + str(extractRegionStart) + ":" + str(extractRegionEnd) + ":" + "0" + ":" + "+"
                ]
            )
        )
        intergenicFlankBedFile.write(extractRegionExtendedBedFileLine + "\n")

        cnt += 1
        if cnt >= numOfsamples:
            raise StopIteration

except StopIteration:
    pass

annotationFile.close()
intergenicFlankBedFile.close()

intergenicFlankFaFileName = workingDir + "intergenicRegions.fa"
command = "bedtools getfasta -name -s -fi " + genome_file + " -bed " + intergenicFlankBedFileName + " -fo " + intergenicFlankFaFileName
print(command)
os.system(command)

for upstreamFlankRegion in flankRegionList:
    cnt = 0

    intergenicFlankFaFileName = workingDir + "intergenicRegions.fa"
    intergenicFlankFaFile = open(intergenicFlankFaFileName, "r")

    intergenicOrfFaFileName = workingDir + "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + ".fa"
    intergenicOrfFaFile = open(intergenicOrfFaFileName, "w")

    intergenicOrfBedFileName = workingDir + "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + ".bed"
    intergenicOrfBedFile = open(intergenicOrfBedFileName, "w")

    try:
        for aLine in intergenicFlankFaFile:
            aLine = aLine.rstrip().upper()

            if aLine[0] == ">":
                header = aLine
                headerCols = header[1:].rstrip().split(":")
                headerName = headerCols[0]
                headerChrom = "chr" + headerCols[1][3:]
                headerStart = int(headerCols[2])
                headerEnd = int(headerCols[3])
                headerScore = headerCols[4]
                headerStrand = headerCols[5]
            else:
                searchPos = 0

                while True:
                    aPos = aLine.find("ATG", searchPos)
                    if aPos == -1:
                        break

                    stopPos = -1
                    for i in range(aPos + 3, len(aLine) - 2, 3):
                        codon = aLine[i:i+3]
                        if codon in stopCodons:
                            stopPos = i
                            break

                    # Αν δεν βρεθεί in-frame stop codon, συνέχισε λίγο μετά το ATG
                    if stopPos == -1:
                        searchPos = aPos + 3
                        continue

                    # Έλεγχος ότι υπάρχει αρκετό upstream
                    if aPos >= upstreamFlankRegion:
                        orfSim = aLine[aPos - upstreamFlankRegion : stopPos + 3]

                        if "N" in orfSim:
                            searchPos = stopPos + 3
                            continue

                        if headerStrand == "+":
                            orfStart = headerStart + aPos - upstreamFlankRegion
                            orfEnd = headerStart + stopPos + 3
                        elif headerStrand == "-":
                            orfStart = headerEnd - (stopPos + 3)
                            orfEnd = headerEnd - aPos + upstreamFlankRegion
                        else:
                            searchPos = stopPos + 3
                            continue

                        intergenicOrfFaFile.write(
                            ">" + headerName + ":" + headerChrom + ":" + str(orfStart) + ":" + str(orfEnd) + ":" + headerScore + ":" + headerStrand + "\n" +
                            orfSim + "\n"
                        )

                        orfBedAnnotation = "\t".join(
                            [headerChrom, str(orfStart), str(orfEnd), headerName, headerScore, headerStrand]
                        )
                        intergenicOrfBedFile.write(orfBedAnnotation + "\n")

                        cnt += 1
                        if cnt >= numOfsamples:
                            raise StopIteration

                    # Το επόμενο ORF θα ψαχτεί μετά το stop του προηγούμενου
                    searchPos = stopPos + 3

    except StopIteration:
        pass

    intergenicFlankFaFile.close()
    intergenicOrfFaFile.close()
    intergenicOrfBedFile.close()