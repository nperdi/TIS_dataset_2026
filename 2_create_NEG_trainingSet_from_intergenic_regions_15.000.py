import sys
sys.path.append('../lib')
#from nMersDB import *
import os
import re
annotationFileName="output/dist_GT_20000_LT_50000/final/6_Human_GRCh38_RefSeq_Curated_distGT_20000_LT_50000_chrFilter_complement_merge_slop_sorted_chrFilter.bed"
workingDir="output/NEGATIVE/"
command="mkdir -p "+workingDir
os.system(command)

def find_all(a_str, sub):
    start = 0
    aList=[]
    while True:
        start = a_str.find(sub, start)
        #print (start)
        if start == -1:
            return aList
        else:
            aList+=[start]
            #print (aList)
        start += len(sub) # use start += 1 to find overlapping matches
        
##list(find_all('spam spam spam spam', 'spam')) # [0, 5, 10, 15]
    
###########################################################

refGenome="hg38"
if refGenome=="hg38":
    genome_file="data/GENOME/hg38.fa"
    chromSizes="data/GENOME/hg38.chrom.sizes"
       
###########################################################
numOfsamples=15000
relativeStart=0
relativeEnd=30000
flankRegionList=[0,100]
downstreamLength=300
###########################################################

annotationFile=open(annotationFileName,"r")

negativeFlankBedFileName=workingDir+"negativeRegions.bed"
negativeFlankBedFile=open(negativeFlankBedFileName,"w")



cnt=0
try:
    for aLine in annotationFile:
        if len(aLine)<=1:
            continue

        cols=aLine.rstrip().split("\t")
        chrom="chr"+cols[0][3:]
        
        #######################
        #if chrom=="chr1":
        #    continue
        #######################
    
        start=int(cols[1])
        end=int(cols[2])        
        name=cols[3]
        if end - start<relativeEnd:
            relativeEnd=end - start
        extractRegionStart=start+relativeStart
        extractRegionEnd=start+relativeEnd   
        extractRegionExtendedBedFileLine="\t".join(map(str,[chrom,extractRegionStart,extractRegionEnd,name+":"+chrom+":"+str(extractRegionStart)+":"+str(extractRegionEnd)+":"+"0"+":"+"+"]))
        negativeFlankBedFile.write(extractRegionExtendedBedFileLine+"\n")
        if cnt>=numOfsamples:
            raise StopIteration
        cnt+=1
     
except StopIteration:
    pass

annotationFile.close()
negativeFlankBedFile.close()



negativeFlankFaFileName=workingDir+"negativeRegions.fa"
command="bedtools getfasta -name -s -fi "+genome_file+" -bed "+negativeFlankBedFileName+" -fo "+negativeFlankFaFileName;
print (command)
os.system(command)

for upstreamFlankRegion in flankRegionList:
    cnt=0
    invalidDict={}
    negativeFlankFaFileName=workingDir+"negativeRegions.fa"
    negativeFlankFaFile=open(negativeFlankFaFileName,"r")
    
    negativeSorfFaFileName=workingDir+"negative_trainingSet_Flank-"+str(upstreamFlankRegion)+".fa"
    negativeSorfFaFile=open(negativeSorfFaFileName,"w")
    
    negativeSorfBedFileName=workingDir+"negative_trainingSet_Flank-"+str(upstreamFlankRegion)+".Bed"
    negativeSorfBedFile=open(negativeSorfBedFileName,"w")
    
    try:
        for aLine in negativeFlankFaFile:
            aLine=aLine.rstrip().upper()
            if  aLine[0]==">":  ##if( $line =~ />(.+)/ ){
                header=aLine
                headerCols=header[1:].rstrip().split(":")
                #print (headerCols)
                headerName=headerCols[0]
                headerChrom="chr"+headerCols[1][3:]
                headerStart=int(headerCols[2])
                headerEnd=int(headerCols[3])
                headerScore=headerCols[4]
                headerStrand=headerCols[5]      
            else:
                startCodonList=find_all(aLine,"ATG")
                #print (aLine)
                #print (startCodonList),
                prevpos=0
                for aPos in startCodonList:
                    if aPos-prevpos<500:
                        continue
                    prevpos=aPos
                    if ((aPos>upstreamFlankRegion) and (aPos<len(aLine)-501)):
                        orfSim=aLine.upper()[aPos-upstreamFlankRegion:aPos+downstreamLength]
                        if  "N" in orfSim or "n" in orfSim :
                            #print (orfSim),
                            continue
                        #print (tisConcensusSeq)
                        if headerStrand=="+":
                            orfStart=headerStart+aPos-upstreamFlankRegion
                            orfEnd=headerStart+aPos+downstreamLength
                        if headerStrand=="-":
                            orfStart=headerEnd-aPos-downstreamLength
                            orfEnd=headerEnd-aPos+upstreamFlankRegion
                                     
                        negativeSorfFaFile.write(">"+headerName+":"+headerChrom+":"+str(orfStart)+":"+str(orfEnd)+":"+headerScore+":"+headerStrand+"\n"+orfSim+"\n")
                        orfBedAnnotation="\t".join([headerChrom,str(orfStart),str(orfEnd),headerName,headerScore,headerStrand])
                        #print (orfBedAnnotation)
                        negativeSorfBedFile.write(orfBedAnnotation+"\n")
                        cnt+=1
                        if cnt>=numOfsamples:
                            raise StopIteration
    except StopIteration:
        pass
    
    negativeFlankFaFile.close()
    negativeSorfFaFile.close()
    negativeSorfBedFile.close()