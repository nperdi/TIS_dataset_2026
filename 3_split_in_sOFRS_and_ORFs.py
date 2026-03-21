import os

inputWorkingDir = "output/intergenic_ORFs/"
flankRegionList = [0, 100]
lengthThreshold = 300

baselineNegativeSorfDir = "output/negative_sORFs_lt" + str(lengthThreshold) + "/"
longIntergenicOrfDir = "output/intergenic_ORFs_ge" + str(lengthThreshold) + "/"

os.makedirs(baselineNegativeSorfDir, exist_ok=True)
os.makedirs(longIntergenicOrfDir, exist_ok=True)

for upstreamFlankRegion in flankRegionList:
    inputOrfBedFileName = inputWorkingDir + "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + ".bed"
    inputOrfFaFileName = inputWorkingDir + "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + ".fa"

    if not os.path.exists(inputOrfBedFileName):
        print("Missing BED:", inputOrfBedFileName)
        continue

    if not os.path.exists(inputOrfFaFileName):
        print("Missing FASTA:", inputOrfFaFileName)
        continue

    negativeSorfBedFileName = os.path.join(
        baselineNegativeSorfDir,
        "baseline_negative_sORFs_Flank-" + str(upstreamFlankRegion) + "_lt" + str(lengthThreshold) + ".bed"
    )
    longIntergenicOrfBedFileName = os.path.join(
        longIntergenicOrfDir,
        "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + "_ge" + str(lengthThreshold) + ".bed"
    )

    negativeSorfFaFileName = os.path.join(
        baselineNegativeSorfDir,
        "baseline_negative_sORFs_Flank-" + str(upstreamFlankRegion) + "_lt" + str(lengthThreshold) + ".fa"
    )
    longIntergenicOrfFaFileName = os.path.join(
        longIntergenicOrfDir,
        "intergenic_ORFs_Flank-" + str(upstreamFlankRegion) + "_ge" + str(lengthThreshold) + ".fa"
    )

    recordCategoryList = []
    negativeSorfBedCount = 0
    longIntergenicOrfBedCount = 0

    ####################################################################
    # 1. Split BED
    ####################################################################
    with open(inputOrfBedFileName, "r") as inputBedFile, \
         open(negativeSorfBedFileName, "w") as negativeSorfBedFile, \
         open(longIntergenicOrfBedFileName, "w") as longIntergenicOrfBedFile:

        for aLine in inputBedFile:
            aLine = aLine.rstrip()
            if len(aLine) <= 1:
                continue

            cols = aLine.split("\t")
            if len(cols) < 3:
                continue

            start = int(cols[1])
            end = int(cols[2])
            orfLength = end - start

            if orfLength < lengthThreshold:
                negativeSorfBedFile.write(aLine + "\n")
                recordCategoryList.append("negative_sorf")
                negativeSorfBedCount += 1
            else:
                longIntergenicOrfBedFile.write(aLine + "\n")
                recordCategoryList.append("long_orf")
                longIntergenicOrfBedCount += 1

    ####################################################################
    # 2. Split FASTA using same order as BED
    ####################################################################
    negativeSorfFaCount = 0
    longIntergenicOrfFaCount = 0
    recordIndex = -1
    currentHeader = None

    with open(inputOrfFaFileName, "r") as inputFaFile, \
         open(negativeSorfFaFileName, "w") as negativeSorfFaFile, \
         open(longIntergenicOrfFaFileName, "w") as longIntergenicOrfFaFile:

        for aLine in inputFaFile:
            aLine = aLine.rstrip()
            if len(aLine) <= 0:
                continue

            if aLine[0] == ">":
                recordIndex += 1
                currentHeader = aLine
            else:
                seq = aLine

                if recordCategoryList[recordIndex] == "negative_sorf":
                    negativeSorfFaFile.write(currentHeader + "\n")
                    negativeSorfFaFile.write(seq + "\n")
                    negativeSorfFaCount += 1
                else:
                    longIntergenicOrfFaFile.write(currentHeader + "\n")
                    longIntergenicOrfFaFile.write(seq + "\n")
                    longIntergenicOrfFaCount += 1

    ####################################################################
    # 3. Print totals
    ####################################################################
    print("\n======================================")
    print("Flank region:", upstreamFlankRegion)
    print("Threshold   :", lengthThreshold)
    print("--------------------------------------")
    print("Negative sORF BED (< " + str(lengthThreshold) + "):", negativeSorfBedCount)
    print("Long ORF BED (>= " + str(lengthThreshold) + "):", longIntergenicOrfBedCount)
    print("Negative sORF FA  (< " + str(lengthThreshold) + "):", negativeSorfFaCount)
    print("Long ORF FA   (>= " + str(lengthThreshold) + "):", longIntergenicOrfFaCount)
    print("--------------------------------------")
    print("Negative sORF BED:", negativeSorfBedFileName)
    print("Long ORF BED     :", longIntergenicOrfBedFileName)
    print("Negative sORF FA :", negativeSorfFaFileName)
    print("Long ORF FA      :", longIntergenicOrfFaFileName)
    print("======================================")