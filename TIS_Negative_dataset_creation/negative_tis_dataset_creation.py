"""
Script purpose:
This script extracts candidate negative TIS examples from intergenic genomic regions.

Main steps:
1. Read a BED file containing filtered intergenic regions.
2. Extract genomic subregions from the reference genome for a fixed window
   starting at each intergenic region.
3. Save these regions in BED format and retrieve their DNA sequences
   using bedtools getfasta.
4. Scan each extracted sequence for ATG codons.
5. Around each ATG, build a fixed TIS-centered sequence window consisting of:
   - an upstream context of variable length
   - the ATG codon
   - a fixed downstream context
6. Keep only candidate windows that:
   - fit completely inside the extracted intergenic sequence
   - do not contain ambiguous base 'N'
   - optionally do not contain an in-frame stop codon in the downstream window
7. Save the accepted candidate TIS windows in:
   - FASTA format
   - BED format
8. Repeat the process for multiple upstream context lengths.

Inputs:
- BED file with filtered intergenic regions
- Reference genome FASTA file

Outputs:
- BED file with extracted intergenic regions
- FASTA file with extracted intergenic region sequences
- FASTA files with candidate intergenic negative TIS windows
- BED files with candidate intergenic negative TIS windows

Notes:
- The script is designed to generate negative TIS examples from intergenic DNA.
- Each candidate sequence is centered on an ATG and includes a configurable
  upstream and downstream context.
- If forbidInFrameStopInDownstreamWindow is enabled, candidate windows with
  an in-frame stop codon shortly after the ATG are excluded.
- The number of extracted regions and resulting candidate TIS examples is
  limited by numOfsamples.
"""

import sys
sys.path.append('../lib')
# from nMersDB import *
import os
genome_file = "../data/GENOME/hg38.fa"
annotationFileName = "../intergenic_region_parsing/output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir = "output/intergenic_TIS/"
command = "mkdir -p " + workingDir
os.system(command)

###########################################################


    

###########################################################
numOfsamples = 15000
relativeStart = 0
relativeEnd = 30000

upstreamContextLengthList = [100, 300]
downstreamContextLength = 500

forbidInFrameStopInDownstreamWindow = True
stopCodons = {"TAA", "TAG", "TGA"}
###########################################################

annotationFile = open(annotationFileName, "r")

intergenicRegionBedFileName = workingDir + "intergenicRegions.bed"
intergenicRegionBedFile = open(intergenicRegionBedFileName, "w")

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
        intergenicRegionBedFile.write(extractRegionExtendedBedFileLine + "\n")

        cnt += 1
        if cnt >= numOfsamples:
            raise StopIteration

except StopIteration:
    pass

annotationFile.close()
intergenicRegionBedFile.close()

intergenicRegionFaFileName = workingDir + "intergenicRegions.fa"
command = "bedtools getfasta -name -s -fi " + genome_file + " -bed " + intergenicRegionBedFileName + " -fo " + intergenicRegionFaFileName
print(command)
os.system(command)

for upstreamContextLength in upstreamContextLengthList:
    cnt = 0

    intergenicRegionFaFile = open(intergenicRegionFaFileName, "r")

    intergenicTisFaFileName = (
        workingDir
        + "intergenic_negative_TIS_up"
        + str(upstreamContextLength)
        + "_down"
        + str(downstreamContextLength)
        + ".fa"
    )
    intergenicTisFaFile = open(intergenicTisFaFileName, "w")

    intergenicTisBedFileName = (
        workingDir
        + "intergenic_negative_TIS_up"
        + str(upstreamContextLength)
        + "_down"
        + str(downstreamContextLength)
        + ".bed"
    )
    intergenicTisBedFile = open(intergenicTisBedFileName, "w")

    try:
        for aLine in intergenicRegionFaFile:
            aLine = aLine.rstrip().upper()

            if len(aLine) <= 0:
                continue

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

                    windowStart = aPos - upstreamContextLength
                    windowEnd = aPos + 3 + downstreamContextLength

                    if windowStart < 0 or windowEnd > len(aLine):
                        searchPos = aPos + 3
                        continue

                    if forbidInFrameStopInDownstreamWindow:
                        hasInFrameStop = False
                        scanEnd = min(aPos + 3 + downstreamContextLength, len(aLine))

                        for i in range(aPos + 3, scanEnd - 2, 3):
                            codon = aLine[i:i+3]
                            if codon in stopCodons:
                                hasInFrameStop = True
                                break

                        if hasInFrameStop:
                            searchPos = aPos + 3
                            continue

                    tisWindowSeq = aLine[windowStart:windowEnd]

                    if "N" in tisWindowSeq:
                        searchPos = aPos + 3
                        continue

                    if headerStrand == "+":
                        tisStart = headerStart + windowStart
                        tisEnd = headerStart + windowEnd
                    elif headerStrand == "-":
                        tisStart = headerEnd - windowEnd
                        tisEnd = headerEnd - windowStart
                    else:
                        searchPos = aPos + 3
                        continue

                    tisName = (
                        headerName
                        + ":"
                        + headerChrom
                        + ":"
                        + str(tisStart)
                        + ":"
                        + str(tisEnd)
                        + ":ATGpos="
                        + str(aPos)
                    )

                    intergenicTisFaFile.write(
                        ">"
                        + tisName
                        + ":"
                        + headerScore
                        + ":"
                        + headerStrand
                        + "\n"
                        + tisWindowSeq
                        + "\n"
                    )

                    tisBedAnnotation = "\t".join(
                        [headerChrom, str(tisStart), str(tisEnd), tisName, headerScore, headerStrand]
                    )
                    intergenicTisBedFile.write(tisBedAnnotation + "\n")

                    cnt += 1
                    if cnt >= numOfsamples:
                        raise StopIteration

                    searchPos = aPos + 3

    except StopIteration:
        pass

    intergenicRegionFaFile.close()
    intergenicTisFaFile.close()
    intergenicTisBedFile.close()

    print("Finished upstream =", upstreamContextLength)
    print("  FASTA:", intergenicTisFaFileName)
    print("  BED  :", intergenicTisBedFileName)
    print("  Count:", cnt)