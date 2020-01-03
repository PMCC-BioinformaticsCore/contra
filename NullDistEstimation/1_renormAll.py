#!/usr/bin/python
# Usage: "python 1_renormAll.py %s %s %s %s" % (contraList, renormPath, renormList, debug)

import math
import os
import sys

if len(sys.argv) < 5:
    print("Error: 1_renormAll.py Not enough arguments")
    sys.exit(0)

contraList = sys.argv[1]
renormPath = sys.argv[2]
renormList = sys.argv[3]
debug = sys.argv[4]

if debug == "T":
    a_and_b = open(
        os.path.join(os.path.split(renormList)[0], "debug/1_abstats.txt"), "w"
    )
    a_and_b.write("SampleName\ta_2\tb_2\t\n")


def renorm_file(path_CONTRA_out, path_out):  # Arguments: two paths to files
    f_CONTRA_out = open(path_CONTRA_out, "r")
    f_out = open(path_out, "w")

    # 1. Process CONTRA input, turn it into nested arrays
    lines = f_CONTRA_out.readlines()
    f_CONTRA_out.close()

    headers = lines[0]
    contra = []
    for line in lines[1:]:
        contra.append(line.split("\t"))
    # For each line in contra:
    # 0 = region ID, 6 = Mean.of.LogRatio, 7 = Adjusted.Mean.of.LogRatio, 14 = tumour.rd, 15 = normal.rd, 16 = tumour.rd.ori, 17 = normal.rd.ori

    # 2. Go through data, calculating x_T, N_T, x_N, N_N
    # as well as creating list of genes to remove (big gains)
    removed_regions = []  # Listed by region ID
    remove_threshold = 2  # Arbitrarily chosen - can be adjusted
    ##### PROBLEM: due to threshold might only take some regions of a particular gene when in reality, want all of them??
    # Does this make a difference?
    x_T = N_T = x_N = N_N = 0.0

    # 4 = OriStCoordinate, 5 = OriEndCoordinate

    for line in contra:
        k = (
            float(line[5]) - float(line[4]) + 1
        )  # Width / number of bases. Inclusive, so add 1.
        if float(line[7]) > remove_threshold:
            removed_regions.append(line[0])
            x_T += float(line[16]) * k
            x_N += float(line[17]) * k
        N_T += float(line[16]) * k
        N_N += float(line[17]) * k
        # Use original (c_s) not adjusted!! (d_s)
    # print removed_regions

    # 3. Calculate a_2, b_2
    n = 1 - float(x_N) / N_N
    t = 1 - float(x_T) / N_T
    a_2 = math.sqrt(n / t)
    b_2 = math.sqrt(t / n)  # or 1/a_2

    # 4. Rescale read depths (and recalculate log ratio??)
    # Add to new columns
    # Need to multiply necessary tumour.rd (d_t) by a_2, normal.rd (d_n) by b_2
    headers = (
        "Gene.Sym\tChr\tOriStCoordinate\tOriEndCoordinate\tMean.of.Log.Ratio.readj\n"
    )  # fix headers
    f_out.write(headers)
    for line in contra:
        temp = [line[2], line[3], line[4], line[5]]
        # calculate new read depths (actually you dont need these)
        # temp.append(str(float(line[14]) * a_2))
        # temp.append(str(float(line[15]) * b_2))
        # Recalculate Mean.of.Log.Ratio
        temp.append(
            str(math.log((float(line[14]) * a_2) / (float(line[15]) * b_2), 2)) + "\n"
        )
        # Write to output
        f_out.write("\t".join(temp))

    # 5. Close file
    f_out.close()

    if debug == "T":
        a_and_b.write(
            "\t".join([os.path.basename(path_CONTRA_out), str(a_2), str(b_2) + "\n"])
        )


# get file paths from list
# note that these are paths, not just names: makes it a bit trickier
dirs = open(contraList, "r")
renormGood = open(renormList, "w")
filepaths = []
bad_filepaths = []
for line in dirs.readlines():
    filepaths.append(line.rstrip("\n\r"))
dirs.close()

middle = "/table/"
suffix = ".CNATable.10rd.10bases.20bins.txt"
outsuffix = "_renorm.txt"

# process...
count = 0
for filepath in filepaths:
    count += 1
    if count % 100 == 0 and debug == "T":
        print(f"...Now processing {str(count)}...")
    filename = os.path.basename(filepath.rstrip("/\\ "))  # want to get file
    cfilepath = os.path.join(filepath + middle, filename + suffix)
    # print cfilepath
    outpath = os.path.join(renormPath, filename + outsuffix)
    if os.path.exists(cfilepath):
        renorm_file(cfilepath, outpath)
        renormGood.write(filename + "\n")
    else:
        bad_filepaths.append(filename)

if debug == "T":
    a_and_b.close()
    print(f"Finished processing out to {renormList}")
    print("Bad paths:")
    print(" ".join(bad_filepaths))
