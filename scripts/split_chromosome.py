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
# Last Updated : 03 Sep 2011 15:00PM


import os


def splitByChromosome(destFolder):

    try:
        os.mkdir(destFolder + "chr/")
    except:
        print("folder exist")

    inputfile = destFolder + "sample.BEDGRAPH"
    outputfile = destFolder + "chr/chr1.txt"
    file = open(inputfile, "r")
    output = open(outputfile, "w")
    check = "1"

    for row in file:
        if row[0] == "#":
            continue

        cols = row.split()
        chr = cols[0].strip("chr")
        if chr != check:
            output.close()
            check = chr
            output = open(destFolder + "chr/chr" + check + ".txt", "w")
        output.write(row)

    output.close()


# JLMod
def splitByChromosome3(infile):
    destFolder = os.path.dirname(infile) + "/"

    try:
        os.mkdir(destFolder + "chr/")
    except:
        print("folder exist")

    # inputfile = destFolder + "sample.BEDGRAPH"
    inputfile = infile
    outputfile = destFolder + "chr/chr1.txt"
    file = open(inputfile, "r")
    output = open(outputfile, "w")
    check = "1"

    for row in file:
        cols = row.split()
        chr = cols[0].strip("chr")
        if chr != check:
            output.close()
            check = chr
            output = open(destFolder + "chr/chr" + check + ".txt", "w")
        output.write(row)

    output.close()


def splitByChromosome2(folder):

    try:
        os.mkdir(folder + "target/")
    except:
        print("folder exist")

    inputfile = folder + "geneRefCoverage.txt"
    outputfile = folder + "target/chr1.txt"
    file = open(inputfile, "r")
    output = open(outputfile, "w")
    check = "1"

    for row in file:
        cols = row.split()
        chr = cols[0].strip("chr")
        if chr != check:
            output.close()
            check = chr
            output = open(folder + "target/chr" + check + ".txt", "w")
        output.write(row)

    output.close()
