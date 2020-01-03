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


def target_breakdown(filename, maxRegionSize, targetRegionSize, out_fname):
    from math import ceil

    a = open(filename)
    output = open(out_fname, "w")

    for x in a:
        x = x.split()
        chr = x[0]
        start = int(x[1])
        end = int(x[2])
        info = ""
        if len(x) > 3:
            info = x[3]

        regionSize = end - start

        if regionSize > maxRegionSize:
            seg = ceil(regionSize / float(targetRegionSize))
            limit = ceil(regionSize / float(seg))

            s_end = start
            for i in range(int(seg)):
                s_start = s_end
                s_end = int(start + (limit * (i + 1)))
                if s_end > end:
                    s_end = end
                output.write(
                    "\t".join([chr, str(s_start), str(s_end), info + "(SPLIT)"]) + "\n"
                )
        else:
            output.write("\t".join([chr, str(start), str(end), info]) + "\n")

    output.close()
