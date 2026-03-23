import sys
sys.path.append('../lib')
# from nMersDB import *
import os

annotationFileName = "output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir = "output/intergenic_ORFs/"
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