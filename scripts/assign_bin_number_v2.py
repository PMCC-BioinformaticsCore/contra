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
# Last Updated : 12 October 2011 16:43PM


def assignBin(binNumber, srcFile, binFile, targetList, minExons):
    # def assignBin(binNumber, srcFile, binFile, minExons):
    from math import log

    sourceFH = open(srcFile)
    src = sourceFH.readlines()
    # binOut = open(binFile, "w")

    minExons = int(minExons)
    count = 0
    logcov_list = []

    # Get the Log Coverage List for the normal sample
    for exon in src:
        exon = exon.split()
        cov1 = float(exon[7])
        cov2 = float(exon[8])
        cov = (cov1 + cov2) / 2
        if cov > 0:
            logcov = log(cov, 2)
        else:
            logcov = 0
        logcov_list.append(logcov)
    # print "Len of logcov_list:", len(logcov_list)

    # Calculate the boundaries of the bins
    minLog = min(logcov_list)
    maxLog = max(logcov_list)
    boundary = (maxLog - minLog) / binNumber
    # print "Min, Max, Boundary, BinNumber: ", minLog, maxLog, boundary, binNumber

    # Split exons to different bins
    bin_list = []
    boundary_dict = {}
    for logcov in logcov_list:
        i = 1
        set_boundary = minLog + boundary
        while logcov > set_boundary:
            i += 1
            set_boundary = minLog + (boundary * i)
        # boundary_dict[i] = set_boundary
        bin_list.append(i)

    for i in range(binNumber + 2):
        boundary_dict[i] = minLog + (boundary * i)

    # Check the number of exons for each bin
    # Merge with the adjacent bins if it is too small
    if len(bin_list) < minExons:
        newMinExons = int(len(bin_list) / 3)
        print(
            f"The total amount of exons ({len(bin_list)}) is smaller than the specified minExons ({minExons})\nSetting minExons to one third of the amount of targeted areas ({newMinExons})"
        )
        minExons = newMinExons

    print("before" + str(bin_list))
    for z in range(1, binNumber + 1):
        element = bin_list.count(z)
        print("Bin", z, "has", element, "elements")
        if element < minExons:
            while bin_list.count(z) != 0:
                if z != binNumber:
                    bin_list[bin_list.index(z)] += 1
                else:
                    bin_list[bin_list.index(z)] -= 1

    print("after" + str(bin_list))
    # Check the number of exons in the last bin
    last_bin_number = sorted(set(bin_list))[-1]

    if len(set(bin_list)) > 1:
        second_last_bin = sorted(set(bin_list))[-2]
    else:
        second_last_bin = last_bin_number

    element = bin_list.count(last_bin_number)
    if element < minExons:
        while bin_list.count(last_bin_number) != 0:
            if last_bin_number != 1:
                # bin_list[bin_list.index(last_bin_number)] -= 1
                bin_list[bin_list.index(last_bin_number)] = second_last_bin
            print("checking through bins")
            print(bin_list)

    final_number_of_bin = len(set(bin_list))

    # List out the binning boundaries in a txt file
    boundary_list = [boundary_dict[x] for x in sorted(set(bin_list))]
    i = 1

    # boundary_file = binFile + str(final_number_of_bin) + ".boundaries.txt"
    boundary_file = binFile + str(binNumber) + ".boundaries.txt"
    boOut = open(boundary_file, "w")
    boOut.write("\t".join([str(0), str(minLog)]) + "\n")
    for z in boundary_list:
        if i == final_number_of_bin:
            z = maxLog
        boOut.write("\t".join([str(i), str(z)]) + "\n")
        i += 1
    boOut.close()

    # Re-sort the bin numbers - get rid of gaps
    curr_z = 1
    bin_number_dict = {}
    for z in sorted(set(bin_list)):
        bin_number_dict[z] = curr_z
        curr_z += 1

    # Append the bin number to the original file
    # binFile = binFile + str(final_number_of_bin)+".txt"
    binFile = binFile + str(binNumber) + ".txt"
    binOut = open(binFile, "w")

    for exons in src:
        exon = exons.split()
        id = int(exon[0])
        gene = exon[1]
        exonNumber = exon[5]

        target = targetList[int(id) - 1]
        if target.id == id:
            chr = target.chr
            oriStart = target.start
            oriEnd = target.end

        else:
            print("ERROR...")

        f_bin_number = str(bin_number_dict[bin_list[count]])
        binOut.write(
            "\t".join([exons.strip("\n"), f_bin_number, chr, oriStart, oriEnd])
        )
        binOut.write("\n")
        count += 1

    binOut.close()
    print(
        "End of assign.bin.number.py with {} exons in {} bins".format(
            count, final_number_of_bin
        )
    )

    return final_number_of_bin
