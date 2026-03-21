import sys
sys.path.append('../lib')
# from nMersDB import *
import os

annotationFileName = "output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir = "output/NEGATIVE/"
command = "mkdir -p " + workingDir
os.system(command)

###########################################################

refGenome = "hg38"
if refGenome == "hg38":
    genome_file = "data/GENOME/hg38.fa"
    chromSizes = "data/GENOME/hg38.chrom.sizes"

###########################################################
numOfsamples = 15000
relativeStart = 0
relativeEnd = 30000
flankRegionList = [0, 100]
stopCodons = {"TAA", "TAG", "TGA"}
###########################################################

annotationFile = open(annotationFileName, "r")

negativeFlankBedFileName = workingDir + "negativeRegions.bed"
negativeFlankBedFile = open(negativeFlankBedFileName, "w")

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
        negativeFlankBedFile.write(extractRegionExtendedBedFileLine + "\n")

        cnt += 1
        if cnt >= numOfsamples:
            raise StopIteration

except StopIteration:
    pass

annotationFile.close()
negativeFlankBedFile.close()

negativeFlankFaFileName = workingDir + "negativeRegions.fa"
command = "bedtools getfasta -name -s -fi " + genome_file + " -bed " + negativeFlankBedFileName + " -fo " + negativeFlankFaFileName
print(command)
os.system(command)

for upstreamFlankRegion in flankRegionList:
    cnt = 0

    negativeFlankFaFileName = workingDir + "negativeRegions.fa"
    negativeFlankFaFile = open(negativeFlankFaFileName, "r")

    negativeSorfFaFileName = workingDir + "negative_trainingSet_Flank-" + str(upstreamFlankRegion) + ".fa"
    negativeSorfFaFile = open(negativeSorfFaFileName, "w")

    negativeSorfBedFileName = workingDir + "negative_trainingSet_Flank-" + str(upstreamFlankRegion) + ".Bed"
    negativeSorfBedFile = open(negativeSorfBedFileName, "w")

    try:
        for aLine in negativeFlankFaFile:
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

                        negativeSorfFaFile.write(
                            ">" + headerName + ":" + headerChrom + ":" + str(orfStart) + ":" + str(orfEnd) + ":" + headerScore + ":" + headerStrand + "\n" +
                            orfSim + "\n"
                        )

                        orfBedAnnotation = "\t".join(
                            [headerChrom, str(orfStart), str(orfEnd), headerName, headerScore, headerStrand]
                        )
                        negativeSorfBedFile.write(orfBedAnnotation + "\n")

                        cnt += 1
                        if cnt >= numOfsamples:
                            raise StopIteration

                    # Το επόμενο ORF θα ψαχτεί μετά το stop του προηγούμενου
                    searchPos = stopPos + 3

    except StopIteration:
        pass

    negativeFlankFaFile.close()
    negativeSorfFaFile.close()
    negativeSorfBedFile.close()