# Emergence and interstate spread of highly pathogenic avian influenza A (H5N1) in dairy cattle

This repository will host scripts and describe how to reproduce the computational analysis in:
*Nguyen, et al. **Emergence and interstate spread of highly pathogenic avian influenza A (H5N1) in dairy cattle**. https://www.biorxiv.org/content/10.1101/2024.05.01.591751v1*

This project is in active development. The paper is a preprint and has not been peer reviewed. 

Links and data are being uploaded as provenance is determined. We cannot publicly upload the genetic sequences from GISAID, but alignment files and xml files include accessions that may be downloaded.

All other data has been deposited at NCBI Genbank with a listing of accessions provided as soon as they have been generated.

If you have specific queries, please submit an issue in the git repo, and we will address it as soon as we are able. Thank you.


### Project structure ###
See the respective README.txt files within the folders for details/instructions.
- [figures](manuscript-figures/): PDF files of manuscript figures.
- [transphylo](transphylo-analysis/): scripts and input data for TransPhylo analysis.
- [treetime](treetime/): Input ML trees and treetime output.
- [reassortment](reassortment-analysis/): Scripts and instructions for reproducing the for reassortment analysis with TreeSort.
- [tmrca](tmrca): BEAUTI xmls with sequence data removed, and annotated MCC trees.
- [discrete-state](discrete-state/): BEAUTI xmls and ouytput used in SpreadD3.
- [wgs-phylogenomics](wgs-phylogenomics/): Scripts and data for WGS B3.13 analysis (phylogenomic tMRCA + spillovers).
- [variant-calling](variant-calling/): Variant calling analyses.

#### TransPhylo Analysis ####
The script ```main.r``` in the tranphylo directory contains the workflow used to infer the transmission on hpai in cattle. Below is the docstring with the usage instructions

```bash
usage: ./tranphylo-analysis/main.r [-h] [-d date] [--g-mean gmean]
                                   [--s-mean smean] [--g-std gstd]
                                   [--s-std sstd] [-i iter]
                                   T

Run transmission analyses with TransPhylo

positional arguments:
  T                     Treefile to analyze

options:
  -h, --help            show this help message and exit
  -d date, --date date  Date last sampled
  --g-mean gmean        Gamma distribution mean of generation time (in years)
  --s-mean smean        Gama distribution mean of sampling time density (in
                        years)
  --g-std gstd          Gamma distribution standard deviation of generation
                        time (in years)
  --s-std sstd          Gamma distribution standard deviation of sampling time
                        density(in years)
  -i iter, --iter iter  Number of MCMC iterations
```

When running the script, ensure the correct selection of parameters for generation time and sampling density. These are measure in years and must be scaled accordingly. For example, if your mean generation time is 5 days, set the mean to (5/365)=0.0137. Use the same logic to set your standard deviation.
