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