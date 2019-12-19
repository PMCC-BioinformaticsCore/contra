#!/usr/bin/env python3

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
# Last Updated : 23 July 2012 16:43PM


import os
import argparse
import sys
import subprocess
import shlex
from multiprocessing import Process, Manager

from scripts.assign_bin_number_v2 import *
from scripts.average_count import *
from scripts.cn_apply_threshold import *
from scripts.convert_gene_coordinate import *
from scripts.convert_targeted_regions import *
from scripts.split_chromosome import *
from scripts.vcf_out import *
from scripts.get_chr_length import *
from scripts.count_libsize import *
from scripts.target_breakdown import *

# VERSION
VERSION = "2.0.10"

# Absolute Path
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))


class Params:
    """
    Class for top-level system parameters
    """

    def __init__(self):
        # command-line option definition
        self.parser = argparse.ArgumentParser()
        # we remove the optionals so the required come first
        self.parser._action_groups.pop()

        # create the group of required params
        self.required = self.parser.add_argument_group("required arguments")

        # create the group of optional
        self.optional = self.parser.add_argument_group("optional arguments")

        self.required.add_argument(
            "-t",
            "--target",
            help="Target region definition file [BED Format]",
            required=True,
        )
        self.required.add_argument(
            "-s",
            "--test",
            help="Alignment file for the test sample [BAM/SAM]",
            required=True,
        )
        self.required.add_argument(
            "-c",
            "--control",
            help="Alignment file for the control sample [BAM/SAM]",
            required=True,
        )
        self.optional.add_argument("-f", "--fasta", help="Reference Genome [FASTA]")
        self.required.add_argument(
            "-o",
            "--outFolder",
            help="the output folder path name to store the output of analysis",
            required=True,
        )
        self.optional.add_argument(
            "--numBin",
            help="Numbers of bins to group regions. User can specify multiple experiments with different number of bins (comma separated) [20]",
            default="20",
        )
        self.optional.add_argument(
            "--minReadDepth",
            help="The threshold for minimum read depth for each bases [10]",
            type=int,
            default=10,
        )
        self.optional.add_argument(
            "--minNBases",
            help="The threshold for minimum number of bases for each target regions [10]",
            type=int,
            default=10,
        )
        self.optional.add_argument(
            "--sam",
            help="If the specified, test and control sample are in SAM [False]",
            action="store_true",
            dest="sam",
            default=False,
        )
        self.optional.add_argument(
            "--bed",
            help="if specified, control will be in BED format [False]",
            action="store_true",
            dest="bedInput",
            default=False,
        )
        self.optional.add_argument(
            "--pval",
            help="The p-value threshold for filtering [0.05]. Applies to Adjusted P-Value.",
            type=float,
            dest="pval",
            default=0.05,
        )
        self.optional.add_argument(
            "--sampleName",
            help="The name to be appended to the front of default output name ['']",
            default="",
        )
        self.optional.add_argument(
            "--nomultimapped",
            help="The option to remove multi-mapped reads [False]",
            action="store_true",
            default=False,
        )
        self.optional.add_argument(
            "-p",
            "--plot",
            help="Plots log-ratio distribution for each bin [False]",
            action="store_true",
            default=False,
        )
        self.optional.add_argument(
            "--minExon",
            help="Minimum number of Exons in one bin (if less than this, bin that contains small number of exons"
            + "will be moved to the adjacent bins) [2000] ",
            type=int,
            default=100,
        )
        self.optional.add_argument(
            "--minControlRdForCall",
            help="Minimum control readdepth for call [5]",
            type=int,
            dest="minControl",
            default=5,
        )

        self.optional.add_argument(
            "--minTestRdForCall",
            help="Minimum test readdepth for call [0]",
            type=int,
            dest="minTest",
            default=0,
        )

        self.optional.add_argument(
            "--minAvgForCall",
            help="Minimum average coverage for call [20]",
            type=int,
            dest="minAvg",
            default=20,
        )

        self.optional.add_argument(
            "--maxRegionSize",
            help="Maximum Region Size in target region (for breaking large region into smaller region. By default, maxRegionSize 0 means no breakdown) [0]",
            default=0,
        )

        self.optional.add_argument(
            "--targetRegionSize",
            help="Target Region Size for breakdown (if maxRegionSize is non zero) [200]",
            type=int,
            default=100,
        )

        self.optional.add_argument(
            "-l",
            "--largeDeletion",
            help="if specified, CONTRA will run large deletion analysis (CBS). User must have DNAcopy R-library installed to run the analysis. [False]",
            action="store_true",
            dest="large",
            default=False,
        )

        self.optional.add_argument(
            "--smallSegment",
            help="CBS segment size for calling large variations [1]",
            type=int,
            default=1,
        )

        self.optional.add_argument(
            "--largeSegment",
            help="CBS segment size for calling large variations [25]",
            type=int,
            default=25,
        )

        self.optional.add_argument(
            "--lrCallStart",
            help="Log ratios start range that will be used to call CNV [-0.3]",
            type=float,
            dest="lrs",
            default=-0.3,
        )

        self.optional.add_argument(
            "--lrCallEnd",
            help="Log ratios end range that will be used to call CNV [0.3]",
            type=float,
            dest="lre",
            default=0.3,
        )

        self.optional.add_argument(
            "--passSize",
            help="Size of exons that passed the p-value threshold compare to the original exon size [0.35]",
            type=float,
            default=0.35,
        )

        ###
        self.optional.add_argument(
            "--removeDups",
            help="if specified, will remove PCR duplicates [False]",
            action="store_true",
            default=False,
        )
        self.optional.add_argument("--version", action="version", version=VERSION)

        # required parameters list
        self.ERRORLIST = []

        # change system parameters based on any command line
        options = self.parser.parse_args()
        print(options)

        if options.target:
            self.TARGET = options.target
        else:
            # self.parser.print_help()
            # self.parser.error("--target not supplied")
            self.ERRORLIST.append("target")

        if options.test:
            if not os.path.isfile(options.test):
                self.ERRORLIST.append("test")
            self.TEST = options.test

        else:
            # self.parser.error("--test not supplied")
            self.ERRORLIST.append("test")

        if options.control:
            if not os.path.isfile(options.test):
                self.ERRORLIST.append("control")
            self.CONTROL = options.control
        else:
            # self.parser.error("--control not supplied")
            self.ERRORLIST.append("control")

        if options.outFolder:
            self.OUTFOLDER = options.outFolder
        else:
            # self.parser.error("--outFolder not supplied")
            self.ERRORLIST.append("outfolder")

        if len(self.ERRORLIST) != 0:
            # self.parser.print_help()
            self.parser.error("Missing required parameters: " + str(self.ERRORLIST))

        if options.numBin:
            binsNumber = options.numBin.split(",")
            try:
                self.NUMBIN = [int(j) for j in binsNumber]
            except:
                self.NUMBIN = [20]

        self.MINREADDEPTH = options.minReadDepth
        self.MINNBASES = options.minNBases
        self.SAM = options.sam
        self.PVAL = options.pval

        if options.sampleName:
            self.SAMPLENAME = options.sampleName
        else:
            self.SAMPLENAME = "No-SampleName"

        self.NOMULTIMAPPED = options.nomultimapped
        self.PLOT = options.plot
        self.BEDINPUT = options.bedInput
        self.MINEXON = options.minExon
        self.MINCONTROL = options.minControl
        self.MINTEST = options.minTest
        self.MINAVG = options.minAvg
        self.MAXREGIONSIZE = options.maxRegionSize
        self.TARGETREGIONSIZE = options.targetRegionSize
        self.LARGE = options.large
        self.SMALLSEGMENT = options.smallSegment
        self.LARGESEGMENT = options.largeSegment
        self.LRE = options.lre
        self.LRS = options.lrs
        self.PASSSIZE = options.passSize
        self.REMOVEDUPS = options.removeDups

        if self.BEDINPUT and self.SAM:
            print("cant parse sam and BED input, setting input to BED")
            self.SAM = False

    def repeat(self):
        # params test
        print("target        :", self.TARGET)
        print("test          :", self.TEST)
        print("control       :", self.CONTROL)
        #        print "fasta        :", self.FASTA
        print("outfolder     :", self.OUTFOLDER)
        print("numBin        :", self.NUMBIN)
        print("minreaddepth  :", self.MINREADDEPTH)
        print("minNBases     :", self.MINNBASES)
        print("sam           :", self.SAM)
        print("pval          :", self.PVAL)
        print("sampleName    :", self.SAMPLENAME)
        print("nomultimapped :", self.NOMULTIMAPPED)
        print("plot          :", self.PLOT)
        print("bedInput      :", self.BEDINPUT)
        print("minExon       :", self.MINEXON)
        print("largeDeletion :", self.LARGE)
        print("removeDups    :", self.REMOVEDUPS)


def checkOutputFolder(outF):
    print("Creating Output Folder :"),

    if outF[len(outF) - 1] == "/":
        outF = outF[: len(outF) - 1]

    try:
        os.makedirs(outF, exist_ok=True)
    except:
        raise IOError("Could not create outputFolder")
    #        print "cannot create folder '{}'" %outF
    #        print "if folder already exist, please specify other folder"
    #        sys.exit(1)

    try:
        os.makedirs(outF + "/table", exist_ok=True)
        os.makedirs(outF + "/plot", exist_ok=True)

        # dont need to create parent folder with makedirs
        os.makedirs(outF + "/buf/ctrData/", exist_ok=True)
        os.makedirs(outF + "/buf/testData/", exist_ok=True)
    except:
        raise IOError("[ERROR: CANNOT CREATE SUBFOLDERS]")

    print(" Done.")

    return outF


# BEDINPUT
def countTotalReads3(params, folder):
    tempFileName = folder + "/temp.txt"
    tempReadFile = open(tempFileName, "w")
    libsize = get_libsize(params.BEDINPUT)
    tempReadFile.write(libsize)
    tempReadFile.close()


def countTotalReads(params, folder):
    if "testData" in folder:
        inF = params.TEST
    else:
        inF = params.CONTROL

    # Get Total ReadCount TODO: use
    getreadcount = os.system(
        "samtools view {} | wc -l > {}/temp.txt".format(inF, folder)
    )


def samToBam(samfile, bamfile):
    args = shlex.split("samtools view -bS {} -o {}".format(samfile, bamfile))
    samtobam = subprocess.call(args)

    return bamfile


def removeMultiMapped(inF, newBAM):
    # Get New BAM Files with mapping quality > 0
    args = shlex.split("samtools view -bq 1 {} -o {}".format(inF, newBAM))
    removeMM = subprocess.call(args)
    print("Multi mapped reads removed.")


def removeDups(inF, newBAM):
    # Remove dups
    args = shlex.split("samtools view -b -F 0x400 {} -o {}".format(inF, newBAM))
    removeDupsCall = subprocess.call(args)
    print("Removed PCR duplicates.")


# BEDINPUT
def convertBamSimple(params, folder, targetList, genomeFile):
    if "testData" in folder:
        inF = params.TEST
        print("Converting TEST Sample... ")
    else:
        inF = params.CONTROL
        print("Converting CONTROL Sample... ")

    # Copy file to working folder
    os.system("cp {} {}".format(inF, folder + "sample.BEDGRAPH"))

    # Split Bedgraph by its chromosomes
    splitByChromosome(folder)

    # Slice the coverage files to only cover the targeted regions
    print("Getting targeted regions DOC...")
    convertGeneCoordinate(targetList, folder)

    # LIBSIZE
    libsize = str(get_libsize(folder + "geneRefCoverage2.txt"))
    tempLibSize = open(folder + "/temp.txt", "w")
    tempLibSize.write(libsize)
    tempLibSize.close()

    print("Targeted regions pre-processing: Done")


def convertBam(params, folder, targetList, genomeFile):
    if "testData" in folder:
        inF = params.TEST
        print("Converting TEST Sample... ")
    else:
        inF = params.CONTROL
        print("Converting CONTROL Sample... ")

    # Convert BAM Files to BEDGRAPH
    bedgraph = folder + "sample.BEDGRAPH"
    args = shlex.split("genomeCoverageBed -ibam {} -bga -g {}".format(inF, genomeFile))

    # output = subprocess.Popen(args, stdout = subprocess.PIPE).communicate()[0]
    iOutFile = open(bedgraph, "w")
    # iOutFile.write(output)
    output = subprocess.Popen(args, stdout=iOutFile).wait()
    iOutFile.close()

    # Split Bedgraph by its chromosomes
    splitByChromosome(folder)

    # Slice the coverage files to only cover the targeted regions
    print("Getting targeted regions DOC...")
    convertGeneCoordinate(targetList, folder)

    # LIBSIZE
    libsize = str(get_libsize(folder + "geneRefCoverage2.txt"))
    tempLibSize = open(folder + "temp.txt", "w")
    tempLibSize.write(libsize)
    tempLibSize.close()

    print("Targeted regions pre-processing: Done")


def analysisPerBin(params, num_bin, outFolder, targetList):
    import shutil

    bufLoc = outFolder + "/buf"
    # Assign bin number to the median and average file
    numBin = assignBin(
        num_bin, bufLoc + "/average.txt", bufLoc + "/bin", targetList, params.MINEXON
    )

    print("Significance Test ...  ")
    rScriptName = os.path.join(scriptPath, "scripts", "cn_analysis.v3.R")
    args = shlex.split(
        "Rscript {} {} {} {} {} {} {} {} {} {} {}".format(
            rScriptName,
            num_bin,
            params.MINREADDEPTH,
            params.MINNBASES,
            outFolder,
            params.SAMPLENAME,
            params.PLOT,
            numBin,
            params.MINCONTROL,
            params.MINTEST,
            params.MINAVG,
        )
    )
    rscr = subprocess.call(args)

    print("Generating Output Files ... ")
    # Analysis of CNV
    tNameList = os.listdir(outFolder + "/table/")
    if num_bin > 1:
        tNameId = str(num_bin) + "bins"
    else:
        tNameId = str(num_bin) + "bin"

    for tName in tNameList:
        if tNameId in tName:
            break

    if "CNATable" in tName:
        tName = tName[: len(tName) - 4]
        tableName = outFolder + "/table/" + tName
        bufTable = bufLoc + "/" + tName
        applyThreshold(
            tableName, bufTable, params.PVAL, 100000
        )  # params.MAXGAP = 100000

        # Large Region CBS
        if params.LARGE != "False":
            print("DEBUG 266a")
            rScriptName2 = os.path.join(scriptPath, "scripts", "large_region_cbs.R")
            args = shlex.split(
                "Rscript {} {} {} {} {} {} {} {} {}".format(
                    rScriptName2,
                    tableName + ".txt",
                    params.SMALLSEGMENT,
                    params.LARGESEGMENT,
                    params.PVAL,
                    params.PASSSIZE,
                    params.LRS,
                    params.LRE,
                    bufLoc,
                )
            )
            rscr2 = subprocess.call(args)
            print(str(args))
        else:
            print("DEBUG 266b")

        # Generate the DNA sequence (for VCF file)
        print("Skipping VCF generation.. use tabular file instead.")

    else:
        print("Table not found")


def removeTempFolder(tempFolderPath):
    import shutil

    shutil.rmtree(tempFolderPath)

    print("Temp Folder Removed")


def main():

    # TODO: use argparse
    # if len(sys.argv) == 2 and sys.argv[1] == "--version":
    #     print(VERSION)
    #     sys.exit(1)

    # option handling
    params = Params()
    params.repeat()

    # output folder handling
    outFolder = checkOutputFolder(params.OUTFOLDER)
    bufLoc = outFolder + "/buf"

    # convert target file
    sorted_target = os.path.join(bufLoc, "target.BED")
    os.system("sort -k1,1 -k2n {} > {}".format(params.TARGET, sorted_target))

    # target breakdown
    if params.MAXREGIONSIZE > 0:
        new_target = os.path.join(bufLoc, "target_breakdown.BED")
        target_breakdown(
            sorted_target, params.MAXREGIONSIZE, params.TARGETREGIONSIZE, new_target
        )
        sorted_target = new_target

    targetList = convertTarget(sorted_target)

    # convert sam to bam if -sam specified
    if params.SAM:
        print("Pre-processing SAM files")

        test_bam = bufLoc + "/test.BAM"
        ctr_bam = bufLoc + "/control.BAM"

        samTest = Process(target=samToBam, args=(params.TEST, test_bam))
        if not params.BEDINPUT:
            samCtr = Process(target=samToBam, args=(params.CONTROL, ctr_bam))

        samTest.start()
        if not params.BEDINPUT:
            samCtr.start()

        samTest.join()
        if not params.BEDINPUT:
            samCtr.join()

        params.TEST = test_bam
        if not params.BEDINPUT:
            params.CONTROL = ctr_bam

    # remove multi mapped reads if --nomultimapped is specified
    if params.NOMULTIMAPPED:
        print("Removing multi-mapped reads")
        test_bam = bufLoc + "/test_reliable.BAM"
        ctr_bam = bufLoc + "/control_reliable.BAM"

        bamTest = Process(target=removeMultiMapped, args=(params.TEST, test_bam))
        if not params.BEDINPUT:
            bamCtr = Process(target=removeMultiMapped, args=(params.CONTROL, ctr_bam))
            bamTest.start()

        if not params.BEDINPUT:
            bamCtr.start()
            bamTest.join()

        if not params.BEDINPUT:
            bamCtr.join()
            params.TEST = test_bam
        if not params.BEDINPUT:
            params.CONTROL = ctr_bam

        ###
    # Remove PCR duplicates if --removeDups specified
    if params.REMOVEDUPS:
        print("Removing reads marked as duplicates (PCR)")

        test_bam = bufLoc + "/test_removedups.BAM"
        ctr_bam = bufLoc + "/control_removedups.BAM"

        bamTest = Process(target=removeDups, args=(params.TEST, test_bam))

        if not params.BEDINPUT:
            bamCtr = Process(target=removeDups, args=(params.CONTROL, ctr_bam))

            bamTest.start()
        if not params.BEDINPUT:
            bamCtr.start()
            bamTest.join()
        if not params.BEDINPUT:
            bamCtr.join()

            params.TEST = test_bam
        if not params.BEDINPUT:
            params.CONTROL = ctr_bam

    # Get Chromosome Length
    genomeFile = bufLoc + "/sample.Genome"
    get_genome(params.TEST, genomeFile)

    # spawn bam converting scripts
    pTest = Process(
        target=convertBam, args=(params, bufLoc + "/testData/", targetList, genomeFile)
    )

    # BEDINPUT
    if not params.BEDINPUT:

        cTest = Process(
            target=convertBam,
            args=(params, bufLoc + "/ctrData/", targetList, genomeFile),
        )
    else:
        cTest = Process(
            target=convertBamSimple,
            args=(params, bufLoc + "/ctrData/", targetList, genomeFile),
        )
    # start the processes
    pTest.start()
    cTest.start()

    # wait for all the processes to finish before continuing
    pTest.join()
    cTest.join()

    # Get the read depth count from temporary folder
    for folder in [bufLoc + "/testData/", bufLoc + "/ctrData/"]:
        fh = open(folder + "temp.txt")
        tmpDP = int(fh.readlines()[0].strip("\n"))

        if "testData" in folder:
            t1 = tmpDP
        else:
            n1 = tmpDP

    print("Test file read depth     = ", t1)
    print("Control file read depth     = ", n1)
    print("Pre-processing Completed. ")

    # Get the Average of the Log Ratio
    print("Getting the Log Ratio ... ")
    testListName = bufLoc + "/testData/geneRefCoverage.txt"
    controlListName = bufLoc + "/ctrData/geneRefCoverage.txt"
    avOut = bufLoc + "/average.txt"
    averageCount(
        testListName,
        controlListName,
        avOut,
        t1,
        n1,
        params.MINREADDEPTH,
        params.MINNBASES,
    )

    # Analysis. [Bin, significance test, large deletion, vcf output]
    print("Binning ... ")
    binProc = []
    for numBin in params.NUMBIN:
        binProc.append(
            Process(target=analysisPerBin, args=(params, numBin, outFolder, targetList))
        )

    for proc in binProc:
        proc.start()

    for proc in binProc:
        proc.join()

    # Removed Temp Folder
    removeTempFolder(bufLoc)


if __name__ == "__main__":
    main()
    print("Done... ")
