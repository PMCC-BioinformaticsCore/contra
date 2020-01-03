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


def get_libsize(bedgraph_file):
    bedgraph = open(bedgraph_file)
    libsize = 0
    for line in bedgraph:
        line = line.split()
        chr = line[0]
        start = int(line[1])
        end = int(line[2])
        cov = float(line[3])
        # cov	= int(line[3])

        libsize += (end - start) * cov

    bedgraph.close()
    return int(libsize)
