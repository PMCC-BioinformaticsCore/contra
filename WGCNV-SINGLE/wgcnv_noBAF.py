#!/usr/bin/python

import sys
import os

USAGE = f"{sys.argv[0]} <sampPath> <mouse|human> [outDir] [outPrefix]"

if len(sys.argv) < 3:
    print(USAGE)
    sys.exit(1)


sampPath = sys.argv[1]
samp = os.path.basename(sampPath)

if sys.argv[2] == "human":
    chrlen_f = "hg19_chrlen.txt"
elif sys.argv[2] == "mouse":
    chrlen_f = "mm9_chrlen.txt"
else:
    print(f"Invalid input - {sys.argv[2]}")

if len(sys.argv) > 3:
    outPath = sys.argv[3]
else:
    outPath = "."

if len(sys.argv) > 4:
    outPrefix = sys.argv[4]
else:
    outPrefix = samp

scriptPath = os.path.realpath(os.path.dirname(sys.argv[0]))
rscript = os.path.join(scriptPath, "wgcnv_withBAFplot.R")


contraF = None
for root, dirs, files in os.walk(sampPath):
    if os.path.basename(root) == "table":
        for f in files:
            if f.endswith("bins.txt"):
                if contraF is not None:
                    print("WARNING: multiple contra output files is detected.")
                    break
                contraF = os.path.join(root, f)

        if contraF:
            break

print(contraF)

variantF = "NA"  # os.path.join(sampPath,samp+"_Ensembl_annotated.tsv")
# if not os.path.isfile(variantF):
# 	raise Exception("Variant file not found %s" % variantF)


if not os.path.exists(outPath):
    os.makedirs(outPath)

outp = os.path.join(outPath, outPrefix)
cmd = f"Rscript {rscript} {contraF} {variantF} {chrlen_f} {outp}"
print(cmd)
os.system(cmd)

cmd2 = f"convert {outp}*byChr*png {outp}_byChr.pdf"
os.system(cmd2)
os.system(f"rm {outp}*byChr*png")
