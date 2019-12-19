# ----------------------------------------------------------------------#
# Copyright (c) 2011, Richard Lupat & Jason Li.
#
# > Source License <
# This file is part of CONTRA.
#
#    CONTRA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CONTRA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CONTRA.  If not, see <http://www.gnu.org/licenses/>.
#
#
# -----------------------------------------------------------------------#
# Last Updated : 12 Oct 2011 11:00AM


def applyThreshold(outputName, bufTable, threshold, maxGap):
    srcFile = outputName + ".txt"
    outFile = bufTable + ".LargeVariations.txt"
    bedFile = bufTable + ".BED"
    fFile = outputName + ".DetailsFILTERED.txt"
    ts = float(threshold)

    # Read and open files
    sourceFH = open(srcFile)
    srcTable = sourceFH.readlines()
    outTable = open(outFile, "w")
    bedOut = open(bedFile, "w")
    filteredTable = open(fFile, "w")

    # header
    outTable.write("Chr \tStartCoordinate \tEndCoordinate \tGenes \tGain.Loss \n")
    filteredTable.write(srcTable[0])

    prevChr = ""
    prevStatus = ""
    prevEnd = -1
    genes = []
    chrList = []

    for exons in srcTable:
        exon = exons.split()
        try:
            adjPVal = float(exon[12])
        except:
            continue

        # print (str(adjPVal)+"  <= "+str(ts))
        if adjPVal <= ts:
            chr = exon[3]
            gene = exon[2]
            status = exon[13]
            start = exon[4]
            end = exon[5]

            # For first row
            if prevEnd == -1:
                gap = 0
            else:
                gap = int(prevEnd) - int(start)

            # Write Filtered Table
            filteredTable.write(exons)

            # Write Bed File
            # bedOut.write(chr.strip("chr") +"\t" +start +"\t"+ end+"\t"+
            #   chr.strip("chr")+":"+start+"-"+end+":"+str(adjPVal)+"\n")
            bedOut.write(
                chr
                + "\t"
                + start
                + "\t"
                + end
                + "\t"
                + chr
                + ":"
                + start
                + "-"
                + end
                + ":"
                + str(adjPVal)
                + "\n"
            )

            if prevChr == "" and prevStatus == "":
                if chr not in chrList:
                    print(chr)
                    chrList.append(chr)
            elif (chr == prevChr) and (status == prevStatus) and (gap < maxGap):
                start = prevStart
            else:
                outTable.write(prevChr + "\t" + prevStart + "\t" + prevEnd + "\t")
                for gsym in genes:
                    outTable.write(gsym + ", ")
                outTable.write("\t" + prevStatus + "\n")
                genes = []

            if gene not in genes:
                genes.append(gene)
            prevChr = chr
            prevStatus = status
            prevStart = start
            prevEnd = end
        elif len(genes) > 0:
            outTable.write(prevChr + "\t" + prevStart + "\t" + prevEnd + "\t")
            for gsym in genes:
                outTable.write(gsym + ", ")
            outTable.write("\t" + prevStatus + "\n")
            prevChr = ""
            prevStatus = ""
            genes = []

    if len(genes) > 0:
        outTable.write(prevChr + "\t" + prevStart + "\t" + prevEnd + "\t")
        for gsym in genes:
            outTable.write(gsym + ", ")
        outTable.write("\t" + prevStatus + "\n")

    filteredTable.close()
    bedOut.close()
    outTable.close()
