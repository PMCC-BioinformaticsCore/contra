#!/usr/bin/env python

# 24/09 wrapper script for estimating the null distribution
# on a per exon and per gene basis from a bunch of CONTRA outputs
# based on wgcnv wrapper
# Requires the following files in the same folder as this script:
# myTSCA.bed.annotated.txt

import sys
import os

USAGE = f"{sys.argv[0]} <CONTRA_out_list.txt> <output_folder> <bed> [debug output T/F] [histograms T/F]"
# CONTRA_out_list = list of CONTRA output folder names / path to, inc name

# If not enough arguments, print how to use the file
debug = "F"
histograms = "F"
if len(sys.argv) < 4:
    print(USAGE)
    sys.exit(1)
if len(sys.argv) > 4:
    debug = sys.argv[4]
if len(sys.argv) > 5:
    histograms = sys.argv[5]

# Grab the sample path and name from command line input
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

# Get contra out list and check it exists
# technically vulnerable since only checking at one point...
contraList = sys.argv[1]
if not os.path.exists(contraList):
    print("Error: couldn't find contra list")
    sys.exit(1)

# Create output folders if needed
# Structure: outPath/renorm, and outpath/histograms if required
outPath = sys.argv[2]
bed = sys.argv[3]
renormPath = os.path.join(outPath, "renorm")
histPath = os.path.join(outPath, "histograms")
debugPath = os.path.join(outPath, "debug")
if not os.path.exists(outPath):
    os.makedirs(outPath)
if not os.path.exists(renormPath):
    os.makedirs(renormPath)
if histograms == "T" and not os.path.exists(histPath):
    os.makedirs(histPath)
if debug == "T" and not os.path.exists(debugPath):
    os.makedirs(debugPath)

# Step 1 - 1_renormAll.py
# Takes in the list of all contra files, renormalises them to renormPath, and creates a new
# list renormList
renormList = os.path.join(outPath, "renorm_list.txt")
renorm_cmd = "python %s/1_renormAll.py %s %s %s %s" % (
    scriptPath,
    contraList,
    renormPath,
    renormList,
    debug,
)
if debug == "T":
    print("Step 1: 1_renormAll.py")
    print(renorm_cmd)
os.system(renorm_cmd)

# Step 2 - 2_getExons.py
# Takes in renormPath and renormList, and creates summary tables exonSummary wgSummary (as well as errors if debug)
exonSummary = os.path.join(outPath, "summary_ex.txt")
wgSummary = os.path.join(outPath, "summary_wg.txt")
getExons_cmd = "python %s/2_getExons.py %s %s %s %s %s %s" % (
    scriptPath,
    renormPath,
    renormList,
    exonSummary,
    wgSummary,
    debug,
    bed,
)
if debug == "T":
    print("Step 2: 2_getExons.py")
    print(getExons_cmd)

os.system(getExons_cmd)


# Step 3 - 3_density.R
# Run this twice - once for wholegene, once for exon
# Requires exonSummary wgSummary + paths out
exonThresh = os.path.join(outPath, "thresh_exons.txt")
wgThresh = os.path.join(outPath, "thresh_wg.txt")
density_wg_cmd = (
    f"Rscript {scriptPath}/3_density.R {wgSummary} {wgThresh} {debug} {histograms}"
)
density_ex_cmd = (
    f"Rscript {scriptPath}/3_density.R {exonSummary} {exonThresh} %{debug} {histograms}"
)

if histograms == "T":
    density_wg_cmd = density_wg_cmd + " " + histPath + "/wg.pdf"
    density_ex_cmd = density_ex_cmd + " " + histPath + "/exon.pdf"

if debug == "T":
    print("Step 3: 3_density.R")
    print(density_wg_cmd)
    print(density_ex_cmd)

os.system(density_wg_cmd)
os.system(density_ex_cmd)
