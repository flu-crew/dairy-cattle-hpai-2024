# Emergence and interstate spread of highly pathogenic avian influenza A (H5N1) in dairy cattle in the USA
[![DOI](https://zenodo.org/badge/795082303.svg)](https://doi.org/10.5281/zenodo.14510300)

This repository will host scripts and describe how to reproduce the computational analysis in:
*Nguyen, et al. **Emergence and interstate spread of highly pathogenic avian influenza A (H5N1) in dairy cattle**. https://www.biorxiv.org/content/10.1101/2024.05.01.591751v1*

The paper is a preprint and has not been peer reviewed. We cannot publicly upload the genetic sequences from GISAID, but alignment files and xml files include accessions that may be downloaded. All other data has been deposited at NCBI Genbank with a listing of accessions provided in the Dataset S4 in the manuscript-supplemental folder.

If you have specific queries, please submit an issue in the git repo, or email Tavis Anderson <tavis.anderson@usda.gov>, and we will address it as soon as we are able. Thank you.

### Project structure ###
See the respective README.txt files within the folders for details/instructions.
- [figures](manuscript-figures/): PDF files of manuscript figures.
- [supplemental](manuscript-supplemental/): PDF files and excel spreadsheets from the supplemental material
- [transphylo](transphylo-analysis/): scripts and input data for TransPhylo analysis.
- [treetime](treetime/): Input ML trees and treetime output.
- [reassortment](reassortment-analysis/): Scripts and instructions for reproducing the reassortment analysis with TreeSort.
- [tmrca](tmrca): BEAUTI xmls with sequence data removed, and annotated MCC trees.
- [discrete-state](discrete-state/): BEAUTI xmls and ouytput used in SpreadD3.
- [wgs-phylogenomics](wgs-phylogenomics/): Scripts and data for WGS B3.13 analysis (phylogenomic tMRCA + spillovers).
- [variant-calling](variant-calling/): Variant calling analyses.
