#!/usr/bin/env python3
import os
import sys

# -----------------------------------------------------------------------------------
# Parts based off 10/09/13 version of get_exons_current.py
# Modified to take command line input
# removed deepcopy cos that was stupid
# -----------------------------------------------------------------------------------
# getExons_cmd = "python 2_getExons.py %s %s %s %s %s %s" % (renormPath, renormList, exonSummary, wgSummary, debug, bed)

if sys.argv < 7:
    print("step 2 error not enough arguments")
    sys.exit(1)

renormPath = sys.argv[1]
renormList = sys.argv[2]
exonSummary = sys.argv[3]
wgSummary = sys.argv[4]
debug = sys.argv[5]
bed = sys.argv[6]


# Global dictionary - {<Gene>: {<Exon> : [List of log ratios]}}
# where Exon = startcoord_endcoord (Should be chromosome_startcoord_endcoord but thats not in renorm)
gdict = {}

# Error file (write out errors here)
if debug == "T":
    error_out = os.path.join(
        os.path.split(exonSummary)[0], "debug/2_getExons_errors.txt"
    )
    errors = open(error_out, "w")


# Create gdict entries based on myTSCA.bed.annotated.txt
# Format: chromosome startcoord endcoord gene
# myTSCA_in = "/home/lij/Current_Routine_Support/TSCA_WGCNV/myTSCA.bed.annotated.txt"
# myTSCA_in = "/home/eyeap/haloplex/halo.bed.annotated.txt" # halo bed?
myTSCA_in = bed
myTSCA = open(myTSCA_in, "r")
for line in myTSCA.readlines():
    entries = line.split("\t")
    entries[3] = entries[3].rstrip("\n")
    if entries[3] not in gdict:
        gdict[entries[3]] = {}
    # gdict[entries[3]]['_'.join(entries[1:3])] = []
    gdict[entries[3]][entries[0] + ":" + "_".join(entries[1:3])] = []
myTSCA.close()

# Helper function
def calculate_average(numblist):
    total = 0.0
    for item in numblist:
        total += float(item)
    return total / float(len(numblist))


# Helper function - list may contain "NA" entries
def calculate_average_withmissing(numblist):
    total = 0.0
    count = 0
    for item in numblist:
        if item != "NA":
            total += float(item)
            count += 1
    if count > 0:
        return total / float(count)
    else:
        return "NA"


# Given a CONTRA output file, generate a list of exons with average log ratio per exon.
def process_file(f_in):
    GENE_IND = 0
    CHR_IND = 1
    START_IND = 2
    END_IND = 3
    LR_IND = 4

    f = open(f_in, "r")
    if debug == "T":
        errors.write("\n\n" + f_in + "\n")
    # 0 - Gene Sym
    # 1 - OriStCoord, 2 - OriEndCoord
    # 3 - Mean.of.LogRatio (renormalised)

    # Exon Data Format:
    # XXX<GeneName>: { <ExonNumber>: [StartCoord, EndCoord, MeanLogRatio] }
    # <GeneName>: { <ExonNumber>: [StartCoord, EndCoord, MeanLogRatio,Chr] }
    exons = {}

    # Temp list of all genes in gdict
    all_genes = gdict.keys()

    # Preprocessing
    contra = []
    for line in f.readlines()[1:]:
        newline = line.split("\t")
        name = newline[GENE_IND].split("(")[0]
        newline[GENE_IND] = name
        newline[-1] = newline[-1].rstrip("\n")
        contra.append(newline)

    f.close()

    currgene = contra[0][GENE_IND]
    temp = [
        contra[0][START_IND],
        contra[0][END_IND],
        [contra[0][LR_IND]],
        contra[0][CHR_IND],
    ]  # start end list of averages
    exons[currgene] = {}
    exon_n = 1
    for line in contra[1:]:
        # If its a new gene...
        if line[GENE_IND] != currgene:
            # Add old to exons dictionary, after calculating average log ratio over the exon.
            temp[2] = calculate_average(temp[2])
            exons[currgene][exon_n] = temp
            currgene = line[GENE_IND]
            temp = [line[START_IND], line[END_IND], [line[LR_IND]], line[CHR_IND]]
            if currgene not in exons:  # if exon data gets broken up...
                exon_n = 1
                exons[currgene] = {}
            else:
                exon_n = len(exons[currgene].keys()) + 1

        # If its a new exon...
        elif line[START_IND] != temp[1]:  # ie coordinates don't match up
            # Add old to exons dictionary, after calculating average log ratio over the exon.
            temp[2] = calculate_average(temp[2])
            exons[currgene][exon_n] = temp
            exon_n += 1
            temp = [line[START_IND], line[END_IND], [line[LR_IND]], line[CHR_IND]]
        else:
            # Part of the current exon. change the current end coord and add average to list of averages.
            temp[2].append(line[LR_IND])
            temp[1] = line[END_IND]
    # Finish off current gene - BUGFIX
    temp[2] = calculate_average(temp[2])
    exons[currgene][exon_n] = temp

    # Add appropriate data to global dictionary gdict
    # If gdict not yet populated (empty):
    for gene in exons.keys():
        # print "gene %s " % gene
        # print "gdict %s " % str(gdict)
        # print "gdict[gene] %s " % str(gdict[gene])
        exon_list = gdict[gene].keys()
        for exon in exons[gene].keys():
            start = int(exons[gene][exon][0])
            end = int(exons[gene][exon][1])
            ch = exons[gene][exon][3].replace("chr", "")
            exon_name = ch + ":" + str(start) + "_" + str(end)
            #            print "DEBUG 33333333 %s  in   %s" % ( exon_name , str(exon_list[0:10]))
            if exon_name in exon_list:
                gdict[gene][exon_name].append(exons[gene][exon][2])
                exon_list.remove(exon_name)
            else:
                # Attempt to match it to something in exon_list: up to 20 bases missing from each side is OK
                matched = False
                for possible_match in exon_list:
                    components = possible_match.split(":")[1].split("_")
                    possible_start = int(components[0])
                    possible_end = int(components[1])
                    start_diff = start - possible_start
                    end_diff = possible_end - end
                    if (
                        start_diff < 21
                        and end_diff < 21
                        and start_diff >= 0
                        and end_diff >= 0
                    ):
                        # Can match these two together
                        if debug == "T":
                            errors.write(
                                "\t"
                                + gene
                                + "_"
                                + possible_match
                                + " is matched with incomplete "
                                + exon_name
                                + "\n"
                            )
                        gdict[gene][possible_match].append(exons[gene][exon][2])
                        exon_list.remove(possible_match)
                        matched = True
                if not matched:
                    # Couldn't find a match for it: output error
                    if debug == "T":
                        errors.write(gene + "_" + exon_name + " could not be found\n")
        # Take care of missing data: remaining exons in exonlist
        for missing_exon in exon_list:
            gdict[gene][missing_exon].append("NA")
            if debug == "T":
                errors.write(gene + "_" + missing_exon + " marked as missing\n")
        all_genes.remove(gene)

    for gene in all_genes:
        # These genes had nothing in the input data. So, add NA to each exon of the gene.
        for exon in gdict[gene]:
            gdict[gene][exon].append("NA")
        if debug == "T":
            errors.write(gene + " not found in input data\n")


# get filenames
dirs = open(renormList, "r")
filenames = []
good_filenames = []
bad_filenames = []
for line in dirs.readlines():
    filenames.append(line.rstrip("\n\r"))
dirs.close()

count = 0
for filename in filenames:
    count += 1
    if count % 100 == 0 and debug == "T":
        print(f"...Now processing {str(count)}...")
    filepath = os.path.join(renormPath, filename + "_renorm.txt")
    if os.path.exists(filepath):
        process_file(filepath)
        good_filenames.append(filename)
    else:
        if debug == "T":
            print(f"Error: Could not find {filename}")
        bad_filenames.append(filename)

# write out gdict...
output_name = exonSummary
out = open(output_name, "w")
# Write headers
out.write("Exon Name\t" + "\t".join(good_filenames) + "\n")  # dont have bad headers!
for gene in sorted(gdict.keys()):
    for exon in sorted(gdict[gene].keys()):
        out.write(gene + "_" + str(exon))
        for ratio in gdict[gene][exon]:
            out.write("\t" + str(ratio))
        out.write("\n")

out.close()

if debug == "T":
    print(f"Finished processing out to {output_name}")


for gene in gdict.keys():
    if ";" in gene:  # If gene needs to be split
        # print 'splitting',gene
        components = gene.split(";")
        for component_gene in components:
            if component_gene not in gdict:
                gdict[component_gene] = {}
        for exon in gdict[gene]:
            gdict[components[0]][exon] = gdict[gene][exon][:]
            gdict[components[1]][exon] = gdict[gene][exon][:]
        del gdict[gene]


wholegene = open(wgSummary, "w")
# write headers
wholegene.write("Gene\t" + "\t".join(good_filenames) + "\n")
for gene in sorted(gdict.keys()):
    wholegene.write(gene)
    for i in range(0, len(good_filenames)):  # For each sample...
        sample_temp = []
        for exon in gdict[gene].keys():  # For each exon (in that gene):
            # print len(gdict[gene][exon])
            # print "i = "+str(i)+" len(gdict[gene][exon])"+ str(len(gdict[gene][exon]))
            sample_temp.append(
                gdict[gene][exon][i]
            )  # collect the log ratio of that exon of that sample (of that gene)
        wholegene.write(
            "\t" + str(calculate_average_withmissing(sample_temp))
        )  # whole gene log ratio = average of log ratios across all exons of that gene
    wholegene.write("\n")
wholegene.close()


if debug == "T":
    print(f"Finished writing whole-gene data to {wgSummary}")
    errors.close()
    print(f"Finished writing errors to {error_out}")
