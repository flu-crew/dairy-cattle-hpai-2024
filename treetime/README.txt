This file contains all of the output of the TreeTime runs as well as the starting trees

Each segment was ran through 4 different TreeTime commands:
1. treetime mugration --tree treefile --states metadata.csv --attribute US-State --confidence --outdir treetime_loc_mugration
2. treetime mugration --tree treefile --states metadata.csv --attribute host --confidence --outdir treetime_host_mugration
3. treetime clock --tree treefile --aln alignment --dates metadata.csv --covariation --outdir treetime_clock
4. treetime --tree treefile --aln alignment --dates metadata.csv --confidence --max-iter 30 --covariation --clock-filter 5 --outdir timetree
    - NOTE: NP and NS segments were ran using --clock-filter=0 because this is a known issue in TreeTime
        - https://github.com/neherlab/treetime/issues/227
