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
# Last Updated : 05 October 2011 16:43PM


import subprocess, shlex


def get_genome(srcFile, genomeOut):
    genome = open(genomeOut, "w")

    args = shlex.split("samtools view -H {}".format(srcFile))
    raw_header = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
    headers = str(raw_header).split("\n")

    for header in headers:
        header = header.split("\t")
        if header[0][1:] != "SQ":
            continue

        genome.write(header[1].strip("SN:") + "\t" + header[2].strip("LN:") + "\n")

    genome.close()
