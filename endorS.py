#!/usr/bin/env python3

"""Script to calculate the Percent on Target (aka Endogenous DNA), clonality and percent of duplicates in a sample from samtools flag stats.
It accepts can accept up to three files: pre-quality, post-quality filtering and post-dedup. We recommend
to use all files but you can also use with a combination of any those samtools flagstats.
"""
import re
import sys
import json
import argparse
import textwrap

parser = argparse.ArgumentParser(prog='endorS.py',
    usage='python %(prog)s [-h] [--version] -r [<samplesfile>.stats] -qF [<samplesfile>.stats] -dedup [<samplesfile>.stats]',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        author:
        Aida Andrades Valtueña (aida.andrades[at]gmail.com)
        
        description:
        %(prog)s calculates Percent on Target (aka Endogenous DNA) from samtools flagstat files and print to screen.
        The Percent on Target reported will be different depending on the combination of samtools flagstat provided.
        This program also calculates clonality (aka Cluster Factor) and percent duplicates when the flagstat file after duplicate removal is provided
        Use --output flag to write results to a file
        '''))
parser.add_argument('--raw', '-r',
                    metavar='<samplefile>.stats', 
                    type=str, nargs='?', 
                    help= 'output of samtools flagstat in a txt file, assumes no quality filtering nor adapter removal performed')
parser.add_argument('--qualityfiltered', '-qf',
                    metavar='<samplefile>.stats', 
                    type=str, nargs='?', 
                    help= 'output of samtools flagstat in a txt file, assumes quality filtering has been performed, must be provided with at least one of the options -r or -dedup')
parser.add_argument('--deduplicated', '-dedup',
                    metavar='<samplefile>.stats', 
                    type=str, nargs='?', 
                    help= 'output of samtools flagstat in a txt file, duplicate removal has been performed')
#parser.add_argument('samtoolsfiles', metavar='<samplefile>.stats', type=str, nargs='+',
#                    help='output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering')
parser.add_argument('-v','--version', action='version', version='%(prog)s 0.4')
parser.add_argument('--output', '-o', nargs='?', help='specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none')
parser.add_argument('--name', '-n', nargs='?', help='specify name for the output file. Default: extracted from the first samtools flagstat file provided')
#parser.add_argument('--dedupflagstats', '-d', type=str, nargs='?')
args = parser.parse_args()

#Check if at least one of the samtools flagstat provided

if ((args.raw is None) and (args.qualityfiltered is None) and (args.deduplicated is None)):
    print("ERROR: no samtools flagstat provided, please provide at least one samtools flagstat files with the flags --raw, --qualityfiltered or -deduplicated.\nRun:\npython endorS.py --help \nfor more information on how to run this script")
    sys.exit()

if ((args.raw is None) and (args.deduplicated is None)):
    print("ERROR: only --qualityfiltered samtools flagstat file provided. No stats can be calculated")
    sys.exit()

try:
    with open(args.raw, 'r') as raw:
        contentsRaw = raw.read()
    #Extract number of total reads
    totalReads = float((re.findall(r'^([0-9]+) \+ [0-9]+ in total',contentsRaw))[0])
    #Extract number of mapped reads pre-quality filtering:
    mappedRaw = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped ',contentsRaw))[0])

    #Calculate Percentage on target raw (aka endogenous raw)
    if totalReads == 0.0:
        endogenousRaw = 0.000000
        print("WARNING: no reads in the fastq input, Percent on Target raw (%) set to 0.000000")
    elif mappedRaw == 0.0:
        endogenousRaw = 0.000000
        print("WARNING: no mapped reads, Percent on Target raw (%) set to 0.000000")
    else:
        endogenousRaw = float("{0:.6f}".format(round((mappedRaw / totalReads * 100), 6)))
except:
    print("No samtools flagstat --raw provided. \nWARNING: none of the Percent on Target stats will be calculated")

try:
    with open(args.qualityfiltered, 'r') as qF:
        contentsqF = qF.read()
    #Extract number of mapped reads post-quality filtering:
    mappedPost = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped',contentsqF))[0])
    #Calculation of Percent on Target modified (aka endogenous DNA post-quality filtering):
    if args.raw is not None:
        if totalReads == 0.0:
            endogenousQF = 0.000000
            print("WARNING: no reads in the fastq input, Percent on Target raw (%) set to 0.000000")
        elif mappedPost == 0.0:
            endogenousQF = 0.000000
            print("WARNING: no mapped reads, Percent on Target modified (%) set to 0.000000")
        else:
            endogenousQF = float("{0:.6f}".format(round(( mappedPost / totalReads * 100),6)))
except:
    print("No samtools flagstat --qualityfiltered provided. \nWARNING: Percent on Target Modified stat will not be calculated")

try:
    with open(args.deduplicated, 'r') as deDup:
        contentsdeDup= deDup.read()
    #Extract number of reads pre dedup
    totalPreDedup = float((re.findall(r'^([0-9]+) \+ [0-9]+ in total',contentsdeDup))[0])
    #Extract number of mapped reads pre-quality filtering:
    mappedDedup = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped ',contentsdeDup))[0])
    #Check if --raw provided and calculate Percent on Target PostDedup
    if args.raw is not None:
        if totalReads == 0.0:
            endogenousDeDup = 0.000000
            print("WARNING: no reads in the fastq input, Percent on Target PostDeDup (%) set to 0.000000")
        elif mappedDedup == 0.0:
            endogenousDeDup = 0.000000
            print("WARNING: no mapped reads, Percent on Target PostDeDup (%) set to 0.000000")
        else:
            endogenousDeDup = float("{0:.6f}".format(round((mappedDedup / totalReads * 100),6)))
    #Calculate clonality (aka cluster factor)
    try:
        clonality=float("{0:.6f}".format(round(( totalPreDedup / mappedDedup),6)))
    except ZeroDivisionError:
        clonality = 0
        print("WARNING: non mapped reads postDeDup and/or preDeDup, check your BAM file. Clonality set to 0.000000")
    #Calculate Percentage of Duplicates
    try:
        percentDuplicates=float("{0:.6f}".format(round((((totalPreDedup - mappedDedup) / totalPreDedup) * 100),6)))
    except ZeroDivisionError:
        percentDuplicates = 0
        print("WARNING: non mapped reads postDeDup and/or preDeDup, check your BAM file. Percent Duplicates set to 0.000000")
except:
    print("No samtools flagstat --deduplicated provided. \nWARNING: Percent on Target PostDedup, Clonality and Percent of Duplicates stats will not be calculated!")

#Setting the name depending on the -name flag:
if args.name is not None:
    name = args.name
else:
    #Set up the name based on the first samtools flagstats:
    if args.raw is not None:
        name = str(((args.raw.rsplit(".",1)[0]).rsplit("/"))[-1])
    elif args.qualityfiltered is not None:
        name = str(((args.qualityfiltered.rsplit(".",1)[0]).rsplit("/"))[-1])
    else:
        name = str(((args.deduplicated.rsplit(".",1)[0]).rsplit("/"))[-1])
#print(name)

#Creating output

if args.raw is not None:
    if args.qualityfiltered is not None:
        #All samtools flagstat provided: Percent Target Raw, Percent Target Modified, Percent Target PostDedup, Clonality, Percent duplicates
        if args.deduplicated is not None:
            jsonOutput={
                "id": "endorSpy",
                "plot_type": "generalstats",
                "pconfig": {
                    "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
                    "percent_on_target_modified": { "max": 100, "min": 0, "title": "Percent on Target modified (%)", "format": '{:,.2f}'},
                    "percent_on_target_postdedup": { "max": 100, "min": 0, "title": "Percent on Target modified (%)", "format": '{:,.2f}'},
                    "clonality": { "max": 100, "min": 0, "title": "Clonality", "format": '{:,.2f}'},
                    "percent_duplicates": { "max": 100, "min": 0, "title": "Percent Duplicates (%)", "format": '{:,.2f}'}
                },
                "data": {
                    name : { 
                        "percent_on_target": endogenousRaw,
                        "percent_on_target_modified": endogenousQF,
                        "percent_on_target_postdedup": endogenousDeDup,
                        "clonality": clonality,
                        "percent_duplicates": percentDuplicates
                        }
                },
            }
            print("Percent on Target raw (%):",endogenousRaw)
            print("Percent on Target modified (%):",endogenousQF)
            print("Percent on Tardet PostDedup (%):", endogenousDeDup)
            print("Clonality:",clonality)
            print("Percent Duplicates (%):", percentDuplicates)
        # Raw + QF: Percent Target Raw, Percent Target Modified
        else:
            jsonOutput={
                "id": "endorSpy",
                "plot_type": "generalstats",
                "pconfig": {
                    "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
                    "percent_on_target_modified": { "max": 100, "min": 0, "title": "Percent on Target modified (%)", "format": '{:,.2f}'}
                },
                "data": {
                    name : { "percent_on_target": endogenousRaw, "percent_on_target_modified": endogenousQF}
                },
            }
            print("Percent on Target raw (%):",endogenousRaw)
            print("Percent on Target modified (%):",endogenousQF)
    # Raw + Dedup: Percent Target Raw, Percent Target PostDedup, Clonality, Percent duplicates
    elif args.deduplicated is not None:
        jsonOutput={
                "id": "endorSpy",
                "plot_type": "generalstats",
                "pconfig": {
                    "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
                    "percent_on_target_postdedup": { "max": 100, "min": 0, "title": "Percent on Target modified (%)", "format": '{:,.2f}'},
                    "clonality": { "max": 100, "min": 0, "title": "Clonality", "format": '{:,.2f}'},
                    "percent_duplicates": { "max": 100, "min": 0, "title": "Percent Duplicates (%)", "format": '{:,.2f}'}
                },
                "data": {
                    name : { 
                        "percent_on_target": endogenousRaw,
                        "percent_on_target_postdedup": endogenousDeDup,
                        "clonality": clonality,
                        "percent_duplicates": percentDuplicates
                        }
                },
            }
        print("Percent on Target raw (%):",endogenousRaw)
        print("Percent on Tardet PostDedup (%):", endogenousDeDup)
        print("Clonality:",clonality)
        print("Percent Duplicates (%):", percentDuplicates)
    
    #Only raw: Percent Target Raw
    else:
        jsonOutput={
            "id": "endorSpy",
            "plot_type": "generalstats",
            "pconfig": 
            {
                "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
            },
            "data": {
                name : { "percent_on_target": endogenousRaw}
            },
        }
        print("Percent on Target raw (%):",endogenousRaw)
#Only Dedup or Dedup + QF provided: clonality and percent duplicates reported
else:
    jsonOutput={
                "id": "endorSpy",
                "plot_type": "generalstats",
                "pconfig": {
                    "clonality": { "max": 100, "min": 0, "title": "Clonality", "format": '{:,.2f}'},
                    "percent_duplicates": { "max": 100, "min": 0, "title": "Percent Duplicates (%)", "format": '{:,.2f}'}
                },
                "data": {
                    name : { "clonality": clonality, "percent_duplicates": percentDuplicates }
                },
    }
    print("Clonality:",clonality)
    print("Percent Duplicates (%):", percentDuplicates)
    



#Checking for print to screen argument:
if args.output is not None:
    #Creating file with the named after the name variable:
    # #Writing the json output:
    fileName = name + "_percent_on_target_mqc.json"
    #print(fileName)
    # with open(fileName, "w+") as outfile:
    json.dump(jsonOutput, outfile)
    print(fileName,"has been generated")
print("All done!")
