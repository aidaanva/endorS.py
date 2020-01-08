# endorS.py
endorS.py calculates endogenous DNA from samtools flagstat files and print to screen

## Usage 
```bash
python endorS.py [-h] [--version] <samplesfile>.stats [<samplesfile>.stats]
```

## Author
Aida Andrades Valtueña (aida.andrades[at]gmail.com)

## Description and options
```bash
positional arguments:
  <samplefile>.stats    output of samtools flagstat in a txt file (at least one required). If two files are supplied,
                        the mapped reads of the second file is divided by the total reads in the first, since it
                        assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --output [OUTPUT], -o [OUTPUT]
                        specify a file format for an output file. Options: <json> for a MultiQC json output. Default:
                        none
```
