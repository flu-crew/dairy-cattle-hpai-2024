### Reassortment analysis on H5N1 2.3.4.4b strains from the US, including the strains from the March-April 2024 outbreak in dairy cattle ###
This analysis uses TreeSort: [https://github.com/flu-crew/TreeSort](!https://github.com/flu-crew/TreeSort).

#### File structure ####
- treesort-analysis.log: a TreeSort-generated reassortment report.
- HA.annotated.treesort.d2.5.tre: HA phylogenetic tree with reassortment annotations on branches (the `rea` tag on the branches can be used to visualize reassortment events in FigTree).
- requirements.txt: a list of Python packages required for this analysis. You can install all by running `pip install -r requirements.txt`.
- EPI_ids.txt: lists the GISAID EPI identifiers for the seqeuences used in this analysis.
- process-wgs.py: a Python script that can be used for processing a combined file of GISAID and cattle-outbreak sequences published in this work, to align, trim, and infer phylogenetic trees for each segment.
- descriptor_H5N1_2344b_HA.csv: a TreeSort descriptor file that provides paths to trees and alignments.
- run-iqtree.py: a script that runs IQ-Tree (iqtree should be installed and be callable as `iqtree`) and reverts the changes that IQ-Tree performs on taxa names (e.g., '|' -> '_').
- treetime-root.py: a script that uses TreeTime to root a given tree based on the molecular clock constraints.

#### Step-by-step reassortment analysis: ####

1. Combine the GISAID sequences (Epi IDs can be found in EPI_ids.txt) with the IRMA sequences from this analysis in a joint file called 'GISAID_NVSL_combined.fasta'

2. Run `python process-wgs.py`. This will align, standardize, trim to cds, and build trees for each gene segment.

3. Root the HA tree by running `python treetime-root.py HA.iqtree.tre HA.final.aln`
4. `mv HA.final.aln.rooted.tre HA.iqtree.rooted.tre`

5. Install TreeSort from [https://github.com/flu-crew/TreeSort](!https://github.com/flu-crew/TreeSort).
6. Run `treesort -i descriptor_H5N1_HA.csv -o HA.annotated.treesort.d2.5.tre --dev 2.5`. TreeSort will identify reassortment events along the HA tree and will give a segment-by-segment report. The tree with reassortment annotations will be saved in `HA.annotated.treesort.d2.5.tre`. The TreeSort report and the annotated tree from this analysis are saved as `treesort-analysis.log` and `HA.annotated.treesort.d2.5.tre` respectively.
