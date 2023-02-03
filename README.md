# endorS.py
endorS.py calculates various endogenous DNA and deduplication metrics  from samtools flagstat files and print to screen.

## Usage
```bash
usage: python endorS.py [-h] [--version] -r [<samplesfile>.stats] -qF [<samplesfile>.stats] -dedup [<samplesfile>.stats]
```

## Author
Aida Andrades Valtueña (aida.andrades[at]gmail.com)

## Description and options
```bash
endorS.py calculates Percent on target (aka Endogenous DNA) from samtools flagstat files and print to screen.
The Percent on target reported will be different depending on the combination of samtools flagstat provided.
This program also calculates clonality (aka Cluster Factor) and percent duplicates when the flagstat file after duplicate removal is provided
Use --output flag to write results to a file

positional arguments:
  <samplefile>.stats    output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in
                        the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering

optional arguments:
  -h, --help            show this help message and exit
  --raw [<samplefile>.stats], -r [<samplefile>.stats]
                        output of samtools flagstat in a txt file, assumes no quality filtering nor duplicate removal performed
  --qualityfiltered [<samplefile>.stats], -q [<samplefile>.stats]
                        output of samtools flagstat in a txt file, assumes some form of quality or length filtering has been performed, must be provided with at least one of the
                        options -r or -dedup
  --deduplicated [<samplefile>.stats], -d [<samplefile>.stats]
                        output of samtools flagstat in a txt file, whereby duplicate removal has been performed on the input reads
  -v, --version         show program's version number and exit
  --output [OUTPUT], -o [OUTPUT]
                        specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none
  --name [NAME], -n [NAME]
                        specify name for the output file. Default: extracted from the first samtools flagstat file provided
  --verbose, -e         increase output verbosity

```
