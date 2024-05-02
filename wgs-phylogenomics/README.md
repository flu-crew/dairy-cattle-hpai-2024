### Phylogenomic analysis of the reassortment-free B3.13 strains ###
Concatenation of B3.13 gene segments and consequent phylogenetic and tMRCA analyses.
Requires IQ-Tree 2, Python 3, and TreeTime.

#### File structure ####
- wgs.B3.13.timetree.covar.relaxed.tre: a time-scaled tree of B3.13 genomes. It can be visualized in FigTree.
- dates.tsv: estimated dates for ancestral nodes along the B3.13 tree with confidence intervals.
- requirements.txt: a list of Python packages required for this analysis. You can install all by running pip `install -r requirements.txt`.
- B3.13.wgs.fasta: Concatenated and aligned genomes of B3.13 strains.
- concatenate-wgs.py: a script to align, trim to cds, and concatenate gene segments.
- iqtree.partitions.nexus: a partition file for IQ-Tree that describes different partitions for each gene segment (this allows for a separate GTR+F+R5 model for each gene segment).
- B3.13.tips.dates.csv: metadata file with dates for each strain in the analysis.

#### Steps for the phylogenomic analysis: ####

1. Build a tree from concatenated genomes using IQ-Tree partitioned model</br>
``iqtree -s B3.13.wgs.fasta -spp iqtree.partitions.nexus -B 1000 -bnni``

2. Perform TreeTime analysis to estimate dates on the ancestral nodes of the B3.13 evolutionary tree.</br>
``treetime --tree iqtree.partitions.nexus.treefile --aln B3.13.wgs.fasta --dates B3.13.tips.dates.csv --covariation --confidence --relax 1.0 0 --max-iter 10 --outdir concat-timetree-cov-filtered/``
