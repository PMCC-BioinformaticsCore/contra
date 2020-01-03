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
# Last Updated : 28 Sep 2011 11:00AM

import sys
import math


def getAverage(list1):
    if len(list1) > 0:
        return float(sum(list1)) / len(list1)

    return 0.0


def getStdDev(list1, avg):
    var = 0.0
    for x in list1:
        var += (avg - x) ** 2

    if (len(list1) - 1) > 0:
        var /= len(list1) - 1

    return math.sqrt(var)


def getMinMax(list1):
    length = len(list1)
    if length != 0:
        min = list1[0]
        max = list1[length - 1]
    else:
        min = 0
        max = 0

    return min, max


def getMedian(list1):
    length = len(list1)
    center = int(length / 2)
    if length == 0:
        median = 0
    elif length % 2 == 0:
        median = (list1[center] + list1[center - 1]) / 2
    else:
        median = list1[center]
    return median


def createDataDict(count, list1, r, offset, id_check, exon_check):
    tDict = {}
    tDictOri = {}

    while count < len(list1):
        t = list1[count].split()
        tId = t[5]
        tExon = t[6]

        if (tId != id_check) or (tExon != exon_check):
            return count, tDict, tDictOri

        tStart = int(t[2])
        tEnd = int(t[3])
        tCov = float(t[4]) / r + offset  # GeoMean Normalisation
        tCovOri = float(t[4]) + offset  # without scaling

        # filling dict
        while tStart < tEnd:
            tDict[tStart] = tCov
            tDictOri[tStart] = tCovOri  # without scaling
            tStart += 1

        count += 1

    return count, tDict, tDictOri


def getFactor(val1, val2):
    r = math.sqrt(val1 * val2)
    r1 = val1 / r
    r2 = val2 / r
    return r1, r2


def averageCount(
    tFile, nFile, averageOut, tReadCount, nReadCount, rd_threshold, minNBases
):
    tFH = open(tFile)
    tList = tFH.readlines()

    nFH = open(nFile)
    nList = nFH.readlines()
    # constant & counter
    OFF = 1
    tCount = 0
    nCount = 0

    # create and open files
    output = open(averageOut, "w")

    # Offset and Ratio for Geometric Mean Normalisation
    r1, r2 = getFactor(tReadCount, nReadCount)
    if rd_threshold > 0:
        # OFF = 0
        OFF = 0.5

    # big loop
    while nCount < len(nList):
        # initialisation, get the chr, geneID, geneName
        init = tList[tCount].split()
        initial = init[5]
        _exon = init[6]
        chr = init[1]
        gene = init[0]
        _start = int(init[2])

        # print("initial: " + initial)
        # print("Check_Init.id: ", check_init[5])
        # print("_Exon: ", _exon)
        # print("Check_Init.exon: ", check_init[6])

        # check if t-gene and n-gene refer to the same gene
        check_init = nList[nCount].split()
        if check_init[5] != initial or check_init[6] != _exon:
            print("Initial: ", initial)
            print("Check_Init.id: ", check_init[5])
            print("_Exon: ", _exon)
            print("Check_Init.exon: ", check_init[6])
            print("Error. Comparing different Gene")
            sys.exit(1)

        # create data dictionary for tumour and normal data (per each regions/ exon)
        tCount, tDict, tDictOri = createDataDict(tCount, tList, r1, OFF, initial, _exon)
        nCount, nDict, nDictOri = createDataDict(nCount, nList, r2, OFF, initial, _exon)
        # check number of bases in the both gene dict
        if len(nDict) != len(tDict):
            print("N:", len(nDict))
            print("T:", len(tDict))
            print("Error. Different length of dict")
            sys.exit(1)

        # compare coverage
        count = _start
        _max = max(nDict.keys())
        ratioList = []
        tumourList = []
        normalList = []
        tumourOriList = []
        normalOriList = []
        while count <= _max:
            # get ratio
            if (nDict[count] < rd_threshold) and (tDict[count] < rd_threshold):
                ratio = 0.0
            else:
                if tDict[count] == 0:
                    tDict[count] = 0.5

                ratio = math.log((float(tDict[count]) / nDict[count]), 2)
                tumourList.append(tDict[count])
                tumourOriList.append(tDictOri[count])
                normalList.append(nDict[count])
                normalOriList.append(nDictOri[count])
                ratioList.append(ratio)

            count += 1

        assert len(ratioList) == len(tumourList), "List do not have the same length"

        ratioLen = len(ratioList)

        # get average
        avg = getAverage(ratioList)
        sd = getStdDev(ratioList, avg)
        tumourAvg = str(round(getAverage(tumourList), 3))
        normalAvg = str(round(getAverage(normalList), 3))
        tumourOriAvg = str(round(getAverage(tumourOriList), 3))
        normalOriAvg = str(round(getAverage(normalOriList), 3))

        # get median
        ratioList.sort()
        min_logratio, max_logratio = getMinMax(ratioList)
        print(ratioList)
        median = getMedian(ratioList)

        # write output
        if ratioLen >= minNBases:
            output.write(initial + "\t" + gene + "\t" + str(ratioLen) + "\t")
            output.write(str(round(avg, 3)) + "\t" + str(count) + "\t" + _exon + "\t")
            output.write(str(round(sd, 3)) + "\t" + tumourAvg + "\t" + normalAvg + "\t")
            output.write(tumourOriAvg + "\t" + normalOriAvg + "\t")
            output.write(
                str(round(median, 3)) + "\t" + str(round(min_logratio, 3)) + "\t"
            )
            output.write(str(round(max_logratio, 3)) + "\n")

    output.close()

    # print "End of averageCount.py with the last target = '%s'" %(initial)
