#!/usr/bin/python

import os
import sys
import fnmatch
import shlex
import subprocess
import shutil
import math
from optparse import OptionParser
from scripts.split_chromosome import splitByChromosome3 as splitByChromosome
from scripts.get_chr_length import *
from scripts.convert_targeted_regions import *
from scripts.convert_gene_coordinate import (
    convertGeneCoordinate2 as convertGeneCoordinate,
)
from multiprocessing import Pool


class Params:
    """
	Class for top-level system parameters:
	"""

    def __init__(self):
        def mult_files(option, opt_str, value, parser):
            args = []
            for arg in parser.rargs:
                if arg[0] != "-":
                    args.append(arg)
                else:
                    del parser.rargs[: len(args)]
                    break
            if getattr(parser.values, option.dest):
                args.extend(getattr(parser.values, option.dest))
            setattr(parser.values, option.dest, args)

        # command-line option definition
        self.parser = OptionParser()
        self.parser.add_option(
            "-t",
            "--target",
            help="Target region definition file [REQUIRED] [BED Format]",
            action="store",
            type="string",
            dest="target",
        )
        self.parser.add_option(
            "-f",
            "--files",
            help="Files to be converted to baselines [REQUIRED] [BAM]",
            action="callback",
            callback=mult_files,
            dest="files",
        )
        self.parser.add_option(
            "-o",
            "--output",
            help="Output folder [REQUIRED]",
            action="store",
            type="string",
            dest="output",
        )
        self.parser.add_option(
            "-c",
            "--trim",
            help="Portion of outliers to be removed before calculating average [Default: 0.2]",
            action="store",
            dest="trim",
            default="0.2",
        )
        self.parser.add_option(
            "-n",
            "--name",
            help="Output baseline name [Default: baseline]",
            action="store",
            dest="name",
            default="baseline",
        )

        # required parameters list:
        self.ERRORLIST = []

        # change system parameters based on any command line arguments
        (options, args) = self.parser.parse_args()
        if options.target:
            self.TARGET = options.target
        else:
            self.ERRORLIST.append("target")

        if options.files:
            self.FILES = options.files
        else:
            self.ERRORLIST.append("files")

        if options.output:
            self.OUTPUT = options.output
        else:
            self.ERRORLIST.append("output")

        if options.trim:
            self.TRIM = float(options.trim)

        if options.name:
            self.NAME = options.name

        if len(self.ERRORLIST) != 0:
            self.parser.print_help()
            self.parser.error("Missing required parameters")

    def repeat(self):
        # params test
        print(f"target	: {self.TARGET}")
        print(f"files	: {self.FILES}")
        print(f"output	: {self.OUTPUT}")
        print(f"trim	: {self.TRIM}")
        print(f"name	: {self.NAME}")


# option handling
params = Params()
# params.repeat()
targetFile = params.TARGET
infiles = params.FILES
output_dir = params.OUTPUT

# Debug
print(" ------ baseline.py ------- ")
print(f"Target: {targetFile}")
for files in infiles:
    print(f"File: {files}")
print(f"Output Directory: {output_dir}")


def make_new_directory(outdir):
    print(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)


print(" ----- creating output directory -----")
make_new_directory(output_dir)
outdir = os.path.join(output_dir, "buf")
make_new_directory(outdir)

targetFile2 = os.path.join(outdir, os.path.basename(targetFile) + ".sorted")
os.system("sort -k1,1 -k2n %s > %s" % (targetFile, targetFile2))


def processInFile(infile):
    infilename = os.path.basename(infile)
    s_outdir = os.path.join(outdir, infilename)
    make_new_directory(s_outdir)
    genomeFile = os.path.join(s_outdir, infilename + ".chrsummary")
    get_genome(infile, genomeFile)

    bedgraph = os.path.join(s_outdir, infilename + ".BEDGRAPH")
    args = shlex.split("genomeCoverageBed -ibam %s -bga -g %s" % (infile, genomeFile))
    iOutFile = open(bedgraph, "w")
    output = subprocess.Popen(args, stdout=iOutFile).wait()
    iOutFile.close()

    targetList = convertTarget(targetFile2)

    splitByChromosome(bedgraph)

    convertGeneCoordinate(targetList, os.path.dirname(bedgraph) + "/")
    bedgraph_tgtonly = bedgraph + ".TARGETONLY"
    bedgraph_tgtonly_avg = bedgraph + ".TARGETONLY.AVERAGE"
    os.rename(
        os.path.join(os.path.dirname(bedgraph), "geneRefCoverage.txt"), bedgraph_tgtonly
    )
    os.rename(
        os.path.join(os.path.dirname(bedgraph), "geneRefCoverage_targetAverage.txt"),
        bedgraph_tgtonly_avg,
    )
    shutil.copy(bedgraph_tgtonly, outdir)
    shutil.copy(bedgraph_tgtonly_avg, outdir)


print("----- Processing Files -----")
pool = Pool(5)
pool.map(processInFile, infiles)


# STEP 2 -> Create Union BedGraph
allfiles_names = [x for x in os.listdir(outdir) if x.endswith("TARGETONLY")]
allfiles_path = [os.path.join(outdir, x) for x in allfiles_names]

args = ["unionBedGraphs", "-header", "-i"] + allfiles_path + ["-names"] + allfiles_names

print(str(args))
fo = os.path.join(outdir, "TARGETONLY.union.txt")
foh = open(fo, "w")

subprocess.Popen(args, stdout=foh).wait()
foh.close()

# STEP 3 -> POOL METHOD
TRIM = params.TRIM  # TRIM 	= 0.2
f = fo
fh = open(f)

fo = os.path.join(output_dir, params.NAME + ".pooled2_TRIM" + str(TRIM) + ".txt")
fo_libsize = fo + ".LIBSIZE"


foh = open(fo, "w")
fo_libsize_h = open(fo_libsize, "w")


fh.readline()

Ns = None
nSamples = None

lineCount = 0
for line in fh:
    lineCount += 1
    line_elems = line.rstrip().split("\t")
    rangeLen = int(line_elems[2]) - int(line_elems[1])
    if not Ns:
        nSamples = len(line_elems) - 3
        Ns = [0 for x in range(nSamples)]
    for i in range(nSamples):
        ind = i + 3
        Ns[i] += int(line_elems[ind]) * rangeLen


fh.close()
fh = open(f)
fh.readline()


lineCount = 0
nSamples_exclude = int(math.floor(TRIM * nSamples))


def gm_mean(xs):
    tmpprod = 1
    p = 1.0 / len(xs)
    for x in xs:
        tmpprod = tmpprod * math.pow(x, p)
    return tmpprod


Ns_gmean = gm_mean(Ns)


def meanstdv(x):
    from math import sqrt

    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    for a in x:
        std = std + (a - mean) ** 2
    std = sqrt(std / float(n - 1))
    return mean, std


libsize = 0
for line in fh:
    lineCount += 1
    line_elems = line.rstrip().split("\t")
    rangeLen = int(line_elems[2]) - int(line_elems[1])
    xs_tmp = [int(x) for x in line_elems[3:]]
    xs = [float(xs_tmp[i]) * Ns_gmean / Ns[i] for i in range(len(xs_tmp))]
    xs.sort()
    xs_trimmed = xs[nSamples_exclude : (nSamples - nSamples_exclude)]
    # trimmed_mean=sum(xs_trimmed)/float(len(xs_trimmed))
    trimmed_mean, std = meanstdv(xs_trimmed)
    libsize += rangeLen * trimmed_mean
    foh.write("\t".join(line_elems[0:3] + [str(trimmed_mean), str(std)]) + "\n")


fo_libsize_h.write(str(int(round(libsize))))

fh.close()
foh.close()
fo_libsize_h.close()


def removeTempFolder(tempFolderPath):
    import shutil

    shutil.rmtree(tempFolderPath)
    print("Temp Folder Removed")


# Removed Temp Folder
# removeTempFolder(outdir)
