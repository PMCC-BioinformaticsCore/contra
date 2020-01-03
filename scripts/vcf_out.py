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
# Last Updated : 09 Apr 2011 11:00AM


def vcf_out(inF, outF):
    import math

    f = file.readlines(open(inF))
    vcf = open(outF, "w")

    # header
    vcf.write("##fileformat=VCFv4.0\n")
    vcf.write("##reference=1000GenomesPilot-NCBI36\n")
    vcf.write(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant"\n'
    )
    vcf.write(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record"\n'
    )
    vcf.write('##ALT=<ID=CNV,Description="Copy number variable region"\n')
    vcf.write("#CHROM \tPOS \tID \tREF \tALT \tQUAL \tFILTER \tINFO\n")

    count = 0
    while count < len(f):
        if count % 2 == 0:
            region = f[count].strip(">\n")
            region = region.split(":")
            chr = region[0]

            adjPVal = float(region[2])
            if adjPVal <= 0:
                adjPVal = 0
            else:
                adjPVal = -10 * math.log(adjPVal, 10)
            adjPVal = str(round(adjPVal, 3))
            region[1] = region[1].split("-")
            start = region[1][0]
            end = region[1][1]
        else:
            ref = f[count].strip("\n")
            vcf.write(chr + "\t" + start + "\t" + "." + "\t" + ref + "\t")
            vcf.write("<CNV>" + "\t" + adjPVal + "\t" + "PASS" + "\t")
            vcf.write("SVTYPE=CNV;END=" + end + "\n")
        count += 1

    vcf.close()
