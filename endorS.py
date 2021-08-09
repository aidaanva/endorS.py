#!/usr/bin/env python3

"""Script to calculate the endogenous DNA in a sample from samtools flag stats.
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
     %(prog)s calculates endogenous DNA from samtools flagstat files and print to screen
     Use --output flag to write results to a file
   '''))
parser.add_argument('samtoolsfiles', metavar='<samplefile>.stats', type=str, nargs='+',
                    help='output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering')
parser.add_argument('-v','--version', action='version', version='%(prog)s 0.4')
parser.add_argument('--output', '-o', nargs='?', help='specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none')
parser.add_argument('--name', '-n', nargs='?', help='specify name for the output file. Default: extracted from the first samtools flagstat file provided')
parser.add_argument('--postdedupflagstats', '-p', type=str, nargs='?')
parser.add_argument('--qualityfiltering', '-q', action='store_true', help='calculates endogenous when no quality filtering has been performed')
parser.add_argument('--cfactor', '-c', action='store_true', help='calculates Reads mapped pre deduplication/Read mapped post deduplication. WARNING= Only calculated when 2 or 3 (using the -p option) samtools flagstats provided')
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
        print("WARNING: no reads in the fastq input, Endogenous DNA raw (%) set to 0.000000")
    elif mappedPre == 0.0:
        endogenousPre = 0.000000
        print("WARNING: no mapped reads, Endogenous DNA raw (%) set to 0.000000")
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
        print("WARNING: no reads in the fastq input, Endogenous DNA raw (%) set to 0.000000")
    elif mappedPost == 0.0:
        endogenousPost = 0.000000
        print("WARNING: no mapped reads, Endogenous DNA modified (%) set to 0.000000")
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

#set reads mapped as well as cluster factor to NA
mappedPostD = "NA"
clusterFactor = "NA"

#Read 3rd samtools flag stats if --postdedupflagstats used
if args.postdedupflagstats is not None:
    try:
        #Open the samtools flag stats post-quality filtering:
        with open(args.postdedupflagstats, 'r') as postdedup:
            contentsPostD = postdedup.read()
        #Extract number of mapped reads post-quality filtering:
        mappedPostD = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped',contentsPostD))[0])
        #Calculation of endogenous DNA as post-dedup reads/total reads filtering:
        if totalReads == 0.0:
            print("WARNING: no reads in the fastq input, Endogenous DNA raw (%) set to 0.000000")
        elif mappedPostD == 0.0:
            endogenousPostD = 0.000000
            print("WARNING: no mapped reads, Endogenous DNA modified post dedup (%) set to 0.000000")
        elif mappedPost == "NA":
            print("WARNING: No post quality filtering samtools flag provided, no Endogenous DNA modified post dedup (%) nor cluster factor calculated")
        else:
            endogenousPostD = float("{0:.6f}".format(round((mappedPostD / totalReads * 100),6)))
    except:
        print("No post deduplication samtools flags provided")


 #Calculate cluster factor if -c provided:
if args.cfactor is not False and mappedPost != "NA":
    if args.postdedupflagstats is None:
        clusterFactor = float("{0:.6f}".format(round((mappedPre / mappedPost),6)))
    else:
        clusterFactor = float("{0:.6f}".format(round((mappedPost / mappedPostD),6)))

    

#Check if -p provided, only case when mappedPostD is not NA
if mappedPostD == "NA":
    #Check if second samtools flagstats provided, only then mappedPost is not NA
    if mappedPost == "NA":
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
        },
        "data": {
            name : { "endogenous_dna": endogenousPre}
        },
        }
        print("Endogenous DNA raw (%):",endogenousPre)
    #Check if cluster factor flag provided (-c)
    elif args.cfactor is not False:
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
            "endogenous_dna_post": { "max": 100, "min": 0, "title": "Endogenous DNA Post (%)", "format": '{:,.2f}'},
            "cluster_factor": { "max": 100, "min": 0, "title": "Cluster Factor", "format": '{:,.2f}'}
        },
        "data": {
            name : { "endogenous_dna": endogenousPre, "endogenous_dna_post": endogenousPost, "cluster_factor": clusterFactor}
        ,
        }
        }
        print("Endogenous DNA raw (%):",endogenousPre)
        print("Endogenous DNA modified (%):",endogenousPost)
        print("Cluster factor:",clusterFactor)
    else:
        #Reporting only pre and post quality filtering
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
            "endogenous_dna_post": { "max": 100, "min": 0, "title": "Endogenous DNA Post (%)", "format": '{:,.2f}'}
        },
        "data": {
            name : { "endogenous_dna": endogenousPre, "endogenous_dna_post": endogenousPost}
        },
        }
        print("Endogenous DNA raw (%):",endogenousPre)
        print("Endogenous DNA modified (%):",endogenousPost)
else:
    if args.cfactor is not False:
        #Creating the json file: pre, post, postdedup (-p), cluster factor (-c) reported
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
            "endogenous_dna_post": { "max": 100, "min": 0, "title": "Endogenous DNA Post (%)", "format": '{:,.2f}'},
            "endogenous_dna_post_dedup": { "max": 100, "min": 0, "title": "Endogenous DNA Post Dedup (%)", "format": '{:,.2f}'},
            "cluster_factor": { "max": 100, "min": 0, "title": "Cluster Factor", "format": '{:,.2f}'}
        },
        "data": {
            name : { "endogenous_dna": endogenousPre, "endogenous_dna_post": endogenousPost, "endogenous_dna_post_dedup": endogenousPostD, "cluster_factor": clusterFactor}
        },
        }
        print("Endogenous DNA raw (%):",endogenousPre)
        print("Endogenous DNA modified (%):",endogenousPost)
        print("Endogenous DNA modfied post deduplication (%):",endogenousPostD)
        print("Cluster factor:",clusterFactor)
    else:
        #Creating the json file: pre, post, postdedup (-p) no -c
        jsonOutput={
        "id": "endorSpy",
        "plot_type": "generalstats",
        "pconfig": {
            "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
            "endogenous_dna_post": { "max": 100, "min": 0, "title": "Endogenous DNA Post (%)", "format": '{:,.2f}'},
            "endogenous_dna_post_dedup": { "max": 100, "min": 0, "title": "Endogenous DNA Post Dedup (%)", "format": '{:,.2f}'},
        },
        "data": {
            name : { "endogenous_dna": endogenousPre, "endogenous_dna_post": endogenousPost, "endogenous_dna_post_dedup": endogenousPostD}
        },
        }
        print("Endogenous DNA raw (%):",endogenousPre)
        print("Endogenous DNA modified (%):",endogenousPost)
        print("Endogenous DNA modfied post deduplication (%):",endogenousPostD)

#Checking for print to screen argument:
if args.output is not None:
   #Creating file with the named after the name variable:
   #Writing the json output:
   fileName = name + "_endogenous_dna_mqc.json"
   #print(fileName)
   with open(fileName, "w+") as outfile:
      json.dump(jsonOutput, outfile)
      print(fileName,"has been generated")
print("All done!")
