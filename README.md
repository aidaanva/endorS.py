# endorS.py
endorS.py calculates endogenous DNA from samtools flagstat files and print to screen. Additionally if -c is provided it can calculate the proportion of unique reads mapped when compared to all the reads mapped. 

## Usage
```bash
usage: endorS.py [-h] [--version] <samplesfile>.stats [<samplesfile>.stats]
```

## Author
Aida Andrades Valtueña (aida.andrades[at]gmail.com)

## Description and options
```bash
  endorS.py calculates endogenous DNA from samtools flagstat files and print to screen
  Use --output flag to write results to a file

positional arguments:
  <samplefile>.stats    output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in
                        the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --output [OUTPUT], -o [OUTPUT]
                        specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none
  --name [NAME], -n [NAME]
                        specify name for the output file. Default: extracted from the first samtools flagstat file provided
  --postdedupflagstats [POSTDEDUPFLAGSTATS], -p [POSTDEDUPFLAGSTATS]
  --cfactor, -c         calculates Reads mapped pre deduplication/Read mapped post deduplication. WARNING= Only calculated when 2 or 3 (using the -p option) samtools flagstats provided

```
