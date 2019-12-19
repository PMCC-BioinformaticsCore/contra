# CONTRA
CONTRA: copy number analysis for targeted resequencing.
Li J, Lupat R, Amarasinghe KC, Thompson ER, Doyle MA, Ryland GL, Tothill RW, Halgamuge SK, Campbell IG, Gorringe KL.

Bioinformatics Core Facility, Peter MacCallum Cancer Centre, VIC 3002, Australia. Jason.Li@petermac.org

Method published 2012 in Bioinformatics doi: [10.1093/bioinformatics/bts146]{https://doi.org/10.1093/bioinformatics/bts146}


previously available through
http://contra-cnv.sourceforge.net/


# Changelog
v2.1.0
Uses python3
Improved performance through multithreading and restructuring of loops
Autodetection of file formats (work in progress)

v2.0.8
Included NDE and WGCNV workflows.
Removed FASTA dependency. Removed VCF support. Removed PDF file.
Updated online documentation.

v2.0.7
Added --version

v2.0.6
Added option to remove Duplicates


v2.0.4
Fixed a bug with "chr" named chromosomes in cn_apply_threshold
Catch errors with bam files do not contain reads in a targeted chromosome.


v2.0.3
Fixed baseline.py scalability to large sample size
Now compatable with R2.14+
