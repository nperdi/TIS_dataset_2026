import sys
sys.path.append('../lib')
# from nMersDB import *
import os

annotationFileName = "output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir = "output/intergenic_TIS/"
command = "mkdir -p " + workingDir
os.system(command)

###########################################################

refGenome = "hg38"
if refGenome == "hg38":
    genome_file = "data/GENOME/hg38.fa"

###########################################################
numOfsamples = 15000
relativeStart = 0
relativeEnd = 30000

upstreamContextLengthList = [100, 300]
downstreamContextLength = 50

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