#!/usr/bin/env python3

"""Script to calculate the Percent on Target (aka Endogenous DNA) in a sample from samtools flag stats.
It accepts can accept up to two files: pre-quality and post-quality filtering. We recommend
to use both files but you can also use the pre-quality filtering.
"""
import re
import sys
import json
import argparse
import textwrap

parser = argparse.ArgumentParser(prog='endorS.py',
   usage='python %(prog)s [-h] [--version] <samplesfile>.stats [<samplesfile>.stats]',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Aida Andrades Valtueña (aida.andrades[at]gmail.com)

   description:
     %(prog)s calculates Percent on Target (aka Endogenous DNA) from samtools flagstat files and print to screen
     Use --output flag to write results to a file
   '''))
parser.add_argument('samtoolsfiles', metavar='<samplefile>.stats', type=str, nargs='+',
                    help='output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering')
parser.add_argument('-v','--version', action='version', version='%(prog)s 0.4')
parser.add_argument('--output', '-o', nargs='?', help='specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none')
parser.add_argument('--name', '-n', nargs='?', help='specify name for the output file. Default: extracted from the first samtools flagstat file provided')
parser.add_argument('--dedupflagstats', '-d', type=str, nargs='?')
args = parser.parse_args()

#Open the samtools flag stats pre-quality filtering:
try:
    with open(args.samtoolsfiles[0], 'r') as pre:
        contentsPre = pre.read()
    #Extract number of total reads
    totalReads = float((re.findall(r'^([0-9]+) \+ [0-9]+ in total',contentsPre))[0])
    #Extract number of mapped reads pre-quality filtering:
    mappedPre = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped ',contentsPre))[0])
    #Calculation of endogenous DNA pre-quality filtering:
    if totalReads == 0.0:
        endogenousPre = 0.000000
        print("WARNING: no reads in the fastq input, Percent on Target raw (%) set to 0.000000")
    elif mappedPre == 0.0:
        endogenousPre = 0.000000
        print("WARNING: no mapped reads, Percent on Target raw (%) set to 0.000000")
    else:
        endogenousPre = float("{0:.6f}".format(round((mappedPre / totalReads * 100), 6)))
except:
    print("Incorrect input, please provide at least a samtools flag stats as input\nRun:\npython endorS.py --help \nfor more information on how to run this script")
    sys.exit()
#Check if the samtools stats post-quality filtering have been provided:
try:
    #Open the samtools flag stats post-quality filtering:
    with open(args.samtoolsfiles[1], 'r') as post:
        contentsPost = post.read()
    #Extract number of mapped reads post-quality filtering:
    mappedPost = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped',contentsPost))[0])
    #Calculation of endogenous DNA post-quality filtering:
    if totalReads == 0.0:
        endogenousPost = 0.000000
        print("WARNING: no reads in the fastq input, Percent on Target raw (%) set to 0.000000")
    elif mappedPost == 0.0:
        endogenousPost = 0.000000
        print("WARNING: no mapped reads, Percent on Target modified (%) set to 0.000000")
    else:
        endogenousPost = float("{0:.6f}".format(round((mappedPost / totalReads * 100),6)))
except:
    print("Only one samtools flagstat file provided")
    #Set the number of reads post-quality filtering to 0 if samtools
    #samtools flag stats not provided:
    mappedPost = "NA"

#Setting the name depending on the -name flag:
if args.name is not None:
    name = args.name
else:
    #Set up the name based on the first samtools flagstats:
    name= str(((args.samtoolsfiles[0].rsplit(".",1)[0]).rsplit("/"))[-1])
#print(name)

#set reads mapped as well as Clonality to NA
mappedPostD = "NA"
clusterFactor = "NA"

#Read 3rd samtools flag stats if --dedupflagstats used
if args.dedupflagstats is not None:
    try:
        #Open the samtools flag stats post-quality filtering:
        with open(args.dedupflagstats, 'r') as postdedup:
            contentsPostD = postdedup.read()
        #Extract number of mapped reads post-quality filtering:
        mappedPostD = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped',contentsPostD))[0])
        #Calculation of endogenous DNA as post-dedup reads/total reads filtering:
        if totalReads == 0.0:
            print("WARNING: no reads in the fastq input, Percent on Target raw (%) set to 0.000000")
        elif mappedPostD == 0.0:
            endogenousPostD = 0.000000
            print("WARNING: no mapped reads, Percent on Target modified post dedup (%) set to 0.000000")
        elif mappedPost == "NA":
            print("WARNING: No post quality filtering samtools flag provided, no Percent on Target modified post dedup (%) nor clonality calculated")
        else:
            endogenousPostD = float("{0:.6f}".format(round((mappedPostD / totalReads * 100),6)))
            clusterFactor = float("{0:.6f}".format(round((mappedPre / mappedPostD),6)))
            percentDuplicates = float("{0:.6f}".format(round((mappedPostD / mappedPost * 100),6)))
    except:
        print("No post deduplication samtools flags provided")


#Check if -p provided, only case when mappedPostD is not NA
if mappedPostD == "NA":
    #Check if second samtools flagstats provided, only then mappedPost is not NA
    if mappedPost == "NA":
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
        },
        "data": {
            name : { "percent_on_target": endogenousPre}
        },
        }
        print("Percent on Target raw (%):",endogenousPre)
    else:
        #Reporting only pre and post quality filtering
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
            "percent_on_target_post": { "max": 100, "min": 0, "title": "Percent on Target Post (%)", "format": '{:,.2f}'}
        },
        "data": {
            name : { "percent_on_target": endogenousPre, "percent_on_target_post": endogenousPost}
        },
        }
        print("Percent on Target raw (%):",endogenousPre)
        print("Percent on Target modified (%):",endogenousPost)
else:
    #Creating the json file: pre, post, postdedup (-d)
    jsonOutput={
    "id": "endorSpy",
    "plot_type": "generalstats",
    "pconfig": {
        "percent_on_target": { "max": 100, "min": 0, "title": "Percent on Target (%)", "format": '{:,.2f}'},
        "percent_on_target_post": { "max": 100, "min": 0, "title": "Percent on Target Post (%)", "format": '{:,.2f}'},
        "percent_on_target_post_dedup": { "max": 100, "min": 0, "title": "Percent on Target Post Dedup (%)", "format": '{:,.2f}'},
        "clonality": { "max": 100, "min": 0, "title": "Clonality", "format": '{:,.2f}'},
        "percent_duplicates": { "max": 100, "min": 0, "title": "Percent Duplicates (%)", "format": '{:,.2f}'}
    },
    "data": {
        name : { "percent_on_target": endogenousPre, 
        "percent_on_target_post": endogenousPost, 
        "percent_on_target_post_dedup": endogenousPostD, 
        "clonality": clusterFactor,
        "percent_duplicates": percentDuplicates}
    },
    }
    print("Percent on Target raw (%):",endogenousPre)
    print("Percent on Target modified (%):",endogenousPost)
    print("Percent on Target modified post deduplication (%):",endogenousPostD)
    print("Clonality:",clusterFactor)
    print("Percent Duplicates (%):",percentDuplicates)

#Checking for print to screen argument:
if args.output is not None:
   #Creating file with the named after the name variable:
   #Writing the json output:
   fileName = name + "_percent_on_target_mqc.json"
   #print(fileName)
   with open(fileName, "w+") as outfile:
      json.dump(jsonOutput, outfile)
      print(fileName,"has been generated")
print("All done!")
