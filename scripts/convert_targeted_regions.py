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
# Last Updated : 28 Mar 2011 11:00AM


class Target:
    """
	Class for target regions
	"""

    population = 0

    def __init__(self):
        self.id = 0
        self.gene = "unknown"
        self.chr = "chr1"
        self.start = 0
        self.end = 0
        self.numberExon = 0
        self.oriStart = 0
        self.oriEnd = 0


def convertTarget(target):
    targets = open(target)

    targetList = []

    count = 0
    for region in targets:
        region = region.split()
        chr = "chr" + region[0].strip("chr")
        start = region[1]
        end = region[2]
        try:
            gene = region[3]
        except:
            gene = "unknown"
        count += 1

        aTarget = Target()
        aTarget.id = count
        aTarget.gene = gene
        aTarget.chr = chr
        aTarget.start = start
        aTarget.end = end
        aTarget.numberExon = 1
        aTarget.oriStart = start
        aTarget.oriEnd = end

        targetList.append(aTarget)

    return targetList
