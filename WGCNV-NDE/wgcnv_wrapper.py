#!/usr/bin/python

# 10/10  wrapper script for WGCNV - server edit
# based on wgcnv_noBAF.py by Jason Li
# Requires the following files in the same folder as this script:
# - 1_renorm.py, 2_process.py, 3_wgsummary.R
# - myTSCA.bed.annotated.txt (used in process)
# - thresh_exon.txt, thresh_wholegene.txt (used in wgsummary)

import sys
import os

USAGE = f"python {sys.argv[0]} <contra_output_basedir> <output> <thresh_wg.txt> <thresh_ex.txt> <bed.annotated.txt> [debug output T/F] [human|mouse]"

# If not enough arguments, print how to use the file
if len(sys.argv) < 5:
    print(USAGE)
    sys.exit(1)
if len(sys.argv) < 7:
    debug = "F"
else:
    debug = sys.argv[6]
if len(sys.argv) < 8:
    species = "human"
else:
    species = sys.argv[7]

bedAnnotatedTxt = sys.argv[5]
# Grab the sample path and name from command line input
sampPath = sys.argv[1]
sampName = os.path.basename(sampPath)
scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))

contraF = None
# Attempt to find the appropriate CONTRA table file
for root, dirs, files in os.walk(sampPath):
    # find the table folder
    if os.path.basename(root) == "table":
        for f in files:
            if f.endswith("bins.txt"):
                if contraF is not None:
                    print("WARNING: multiple contra output files detected.")
                    break
                contraF = os.path.join(root, f)
        if contraF:
            break

# If there isnt a CONTRA file, exit
if contraF == None:
    print("ERROR: Could not find CONTRA output file.")
    sys.exit(1)

# Create output directory if needed
outPath = sys.argv[2]
if not os.path.exists(outPath):
    os.makedirs(outPath)


# 1. RENORM.PY
# makes _renorm.txt
# This may not be necessary, depending if renorm already exists... <- i guess just remake it to ensure no tampering?
# Perhaps come back and tidy this up later?
renorm_out = os.path.normpath(os.path.join(outPath, sampName + "_renorm.txt"))
renorm_cmd = f"python {scriptPath}/1_renorm.py {contraF} {renorm_out} {debug}"
if debug == "T":
    print("............")
    print("Step 1: 1_renorm.py")
    print(renorm_cmd)
os.system(renorm_cmd)

# 2. PROCESS.PY
# makes _exon.txt, _wg.txt
process_exon_out = os.path.normpath(os.path.join(outPath, sampName + "_exon.txt"))
process_wg_out = os.path.normpath(os.path.join(outPath, sampName + "_wg.txt"))
process_cmd = f"python {scriptPath}/2_process.py {renorm_out} {process_wg_out} {process_exon_out} {bedAnnotatedTxt} {debug}"
if debug == "T":
    print("............")
    print("Step 2: 2_process.py")
    print(process_cmd)
os.system(process_cmd)

# 3. WGSUMMARY.R
# makes summary table _wgSummary.txt
summary_out = os.path.normpath(os.path.join(outPath, sampName + "_wgSummary.txt"))
ex_summary_out = os.path.normpath(os.path.join(outPath, sampName + "_exSummary.txt"))
thresh_wg = sys.argv[3]
thresh_exon = sys.argv[4]
summary_cmd = f"Rscript {scriptPath}/3_wgsummary.R {scriptPath} {process_wg_out} {process_exon_out} {summary_out} {ex_summary_out} {thresh_wg} {thresh_exon} {debug}"

if debug == "T":
    print("............")
    print("Step 3: 3_wgsummary.R")
    print(summary_cmd)
os.system(summary_cmd)

print(f"Finished: created {summary_out} {ex_summary_out}")


# 4. PLOT.R

outf = os.path.normpath(os.path.join(outPath, sampName + "_plot.pdf"))
if species == "human":
    species_f = os.path.join(scriptPath, "hg19_chrlen.txt")
else:
    raise Exception("species not implemented: " + species)

plotScript = os.path.join(scriptPath, "4_plot.R")
plot_cmd = f"Rscript {plotScript} {ex_summary_out} {summary_out} {sampName} {species} {species_f} {outf}"


print("............")
print("Step 4: 4_plot.R")
print(plot_cmd)
os.system(plot_cmd)
print("Finished: created {outf}")
